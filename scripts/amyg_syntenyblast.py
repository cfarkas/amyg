#!/usr/bin/env python3

import os
import sys
import argparse
import subprocess
import shutil
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO
from tqdm import tqdm
from intervaltree import IntervalTree

##############################################################################
# CONSTANTS (internal defaults)
##############################################################################

IDENTITY_THRESHOLD = 50.0
COVERAGE_THRESHOLD = 50.0
TOTAL_CONTIG_LENGTH = 20000  # For the 50% ratio check
BLAST_DB_NAME = "blast_db"
BLAST_RESULTS_BASENAME = "blast_results.txt"
DETAILED_CSV_BASENAME = "processed_blast_results.csv"
SYNTENY_CSV_BASENAME = "synteny_blocks.csv"
TEMP_BLAST_DIRNAME = "temp_blast_results"
PLOT_PREFIX = "synteny_plots"

##############################################################################
# HELPER FUNCTIONS
##############################################################################

def log(msg):
    print(f"[INFO] {msg}", flush=True)

def ensure_directory_exists(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)

def split_fasta(fasta_file, output_dir, chunk_size):
    """
    Splits the FASTA file into chunks of 'chunk_size' bases each,
    writes them to output_dir as separate .fasta files.
    """
    ensure_directory_exists(output_dir)
    records = SeqIO.parse(fasta_file, "fasta")
    total = sum(1 for _ in SeqIO.parse(fasta_file, "fasta"))  # For tqdm
    pbar = tqdm(total=total, desc="Splitting sequences")
    for record in records:
        seq_length = len(record.seq)
        for start in range(0, seq_length, chunk_size):
            end = min(start + chunk_size, seq_length)
            chunk_seq = record.seq[start:end]
            chunk_id = f"{record.id}_{start}_{end}"
            chunk_record = SeqIO.SeqRecord(chunk_seq, id=chunk_id, description="")
            chunk_path = os.path.join(output_dir, f"{chunk_id}.fasta")
            SeqIO.write(chunk_record, chunk_path, "fasta")
        pbar.update(1)
    pbar.close()

def create_blast_db(chunks_dir, db_path):
    """
    Concatenate all .fasta chunks in 'chunks_dir' into 'all_chunks.fasta',
    then calls makeblastdb to create a DB at 'db_path'.
    """
    input_fasta = os.path.join(chunks_dir, "all_chunks.fasta")
    with open(input_fasta, 'w') as outfile:
        for filename in os.listdir(chunks_dir):
            if filename.endswith('.fasta'):
                with open(os.path.join(chunks_dir, filename)) as infile:
                    outfile.write(infile.read())

    command = f"makeblastdb -in {input_fasta} -dbtype nucl -out {db_path}"
    subprocess.run(command, shell=True, check=True)

def run_blast(chunks_dir, db_path, temp_blast_dir, threads):
    """
    Runs BLASTN for each .fasta chunk in 'chunks_dir' vs. 'db_path'.
    Writes intermediate results to 'temp_blast_dir/xxx_temp.txt',
    merges them into 'temp_blast_dir/concatenated_blast.txt'.
    
    Filters:
      - alignment ratio >= 0.50
      - not self
      - query length >= 0.50 of TOTAL_CONTIG_LENGTH
      - IDENTITY_THRESHOLD, COVERAGE_THRESHOLD
    """
    ensure_directory_exists(temp_blast_dir)
    log(f"Using temporary directory for BLAST results: {temp_blast_dir}")
    files = [f for f in os.listdir(chunks_dir) if f.endswith('.fasta')]

    all_blast_lines = []
    pbar = tqdm(total=len(files), desc="Running BLAST")

    for filename in files:
        query_file = os.path.join(chunks_dir, filename)
        temp_output = os.path.join(temp_blast_dir, f"{filename}_temp.txt")

        # Format 6: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen
        command = (
            f"blastn -query {query_file} -db {db_path} "
            f"-outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend "
            f"sstart send evalue bitscore qlen slen' "
            f"-out {temp_output} "
            f"-num_threads {threads}"
        )
        subprocess.run(command, shell=True, check=True)

        # Filter lines
        with open(temp_output, 'r') as temp_results:
            for line in temp_results:
                parts = line.strip().split()
                qseqid = parts[0]
                sseqid = parts[1]
                identity = float(parts[2])
                align_len = int(parts[3])
                qstart = int(parts[6])
                qend = int(parts[7])
                sstart = int(parts[8])
                send = int(parts[9])
                qlen = int(parts[12])

                alignment_ratio = align_len / qlen
                query_contig_ratio = qlen / TOTAL_CONTIG_LENGTH

                if alignment_ratio >= 0.50 and qseqid != sseqid and query_contig_ratio >= 0.50:
                    coverage = (align_len / qlen) * 100
                    if identity >= IDENTITY_THRESHOLD and coverage >= COVERAGE_THRESHOLD:
                        all_blast_lines.append(line.strip())

        os.remove(temp_output)  # Clean up
        pbar.update(1)
    pbar.close()

    # Write final concatenated results
    concatenated_path = os.path.join(temp_blast_dir, "concatenated_blast.txt")
    with open(concatenated_path, 'w') as outfile:
        for line in all_blast_lines:
            outfile.write(line + "\n")

    log(f"Final concatenated blast output: {concatenated_path}")

def load_and_process_blast_data(blast_file, out_csv):
    """
    Loads the final BLAST file, renames columns, sorts, and writes out_csv.
    Returns a DataFrame.
    """
    cols = [
        'query_seq_id','subject_seq_id','identity_percentage',
        'alignment_length_bp','mismatches_count','gap_opens_count',
        'query_start','query_end','subject_start','subject_end',
        'e_value','bit_score','query_length','subject_length'
    ]
    df = pd.read_csv(blast_file, sep='\t', header=None, names=cols)
    df.drop_duplicates(inplace=True)
    df['query_start'] = df['query_start'].astype(int)
    df.sort_values(by=['query_seq_id','query_start'], inplace=True)
    df.to_csv(out_csv, index=False)
    log(f"Processed BLAST results saved to: {out_csv}")
    return df

def parse_seq_id(seq_id):
    """Extract (contig_name, start, end) from 'contigName_start_end' pattern."""
    match = re.match(r'(.+?)_(\d+)_(\d+)', seq_id)
    if match:
        return match.group(1), int(match.group(2)), int(match.group(3))
    return seq_id, None, None

def process_synteny_blocks(df):
    """
    Merges contiguous regions, removes duplicates, adds columns like duplication_type,
    and returns the final synteny-block DataFrame.
    """
    df[['query_contig','query_seq_start','query_seq_end']] = df['query_seq_id'].apply(
        lambda x: pd.Series(parse_seq_id(x))
    )
    df[['subject_contig','subject_seq_start','subject_seq_end']] = df['subject_seq_id'].apply(
        lambda x: pd.Series(parse_seq_id(x))
    )

    df.dropna(subset=['query_seq_start','subject_seq_start'], inplace=True)
    for c in ['query_seq_start','query_seq_end','subject_seq_start','subject_seq_end']:
        df[c] = df[c].astype(int)

    df.sort_values(by=['query_contig','query_seq_start','subject_contig','subject_seq_start','query_start'], inplace=True)

    merged_rows = []
    for (qc, sc), group in df.groupby(['query_contig','subject_contig']):
        group = group.sort_values(by=['query_seq_start','query_start'])
        current_block = None
        for _, row in group.iterrows():
            if current_block is None:
                current_block = row.copy()
                current_block['identity_percentages_list'] = [row['identity_percentage']]
            else:
                query_adj = (row['query_seq_start'] == current_block['query_seq_end'])
                subj_adj = (row['subject_seq_start'] == current_block['subject_seq_end'])
                if query_adj and subj_adj:
                    current_block['query_seq_end'] = row['query_seq_end']
                    current_block['subject_seq_end'] = row['subject_seq_end']
                    current_block['identity_percentages_list'].append(row['identity_percentage'])
                    current_block['alignment_length_bp'] += row['alignment_length_bp']
                else:
                    current_block['average_percent_identity'] = np.mean(current_block['identity_percentages_list'])
                    merged_rows.append(current_block)
                    current_block = row.copy()
                    current_block['identity_percentages_list'] = [row['identity_percentage']]

        if current_block is not None:
            current_block['average_percent_identity'] = np.mean(current_block['identity_percentages_list'])
            merged_rows.append(current_block)

    merged_df = pd.DataFrame(merged_rows)

    def create_synteny_block_id(r):
        ids = sorted([
            f"{r['query_contig']}_{r['query_seq_start']}_{r['query_seq_end']}",
            f"{r['subject_contig']}_{r['subject_seq_start']}_{r['subject_seq_end']}"
        ])
        return "_vs_".join(ids)

    merged_df['synteny_block_id'] = merged_df.apply(create_synteny_block_id, axis=1)
    merged_df.sort_values(by='alignment_length_bp', ascending=False, inplace=True)
    merged_df.drop_duplicates(subset='synteny_block_id', keep='first', inplace=True)
    merged_df.reset_index(drop=True, inplace=True)

    merged_df['merged_query_seq_id'] = merged_df.apply(
        lambda x: f"{x['query_contig']}_{x['query_seq_start']}_{x['query_seq_end']}", axis=1)
    merged_df['merged_subject_seq_id'] = merged_df.apply(
        lambda x: f"{x['subject_contig']}_{x['subject_seq_start']}_{x['subject_seq_end']}", axis=1)

    merged_df['query_length'] = merged_df['query_seq_end'] - merged_df['query_seq_start']
    merged_df['subject_length'] = merged_df['subject_seq_end'] - merged_df['subject_seq_start']
    merged_df['average_percent_identity'] = merged_df['average_percent_identity'].round(3)

    # duplication_type
    merged_df['duplication_type'] = merged_df['average_percent_identity'].apply(
        lambda x: 'recent' if x >= 95 else 'ancient'
    )

    cols_to_remove = [
        'query_start','query_end','subject_start','subject_end',
        'mismatches_count','gap_opens_count','e_value','bit_score',
        'identity_percentage','identity_percentages_list','synteny_block_id'
    ]
    keep = [c for c in merged_df.columns if c not in cols_to_remove]
    merged_df = merged_df[keep]

    final_cols = [
        'merged_query_seq_id','merged_subject_seq_id',
        'query_contig','query_seq_start','query_seq_end',
        'subject_contig','subject_seq_start','subject_seq_end',
        'average_percent_identity','alignment_length_bp',
        'query_length','subject_length','duplication_type'
    ]
    merged_df = merged_df[final_cols]

    return merged_df

def generate_plots(merged_df, out_prefix):
    merged_df_sorted = merged_df.sort_values(by='average_percent_identity')

    # 1) histogram
    plt.figure(figsize=(10,6))
    plt.hist(merged_df_sorted['average_percent_identity'], bins=20, color='blue', alpha=0.7)
    plt.xlabel('Average Percent Identity')
    plt.ylabel('Frequency')
    plt.title('Distribution of Average Percent Identity')
    plt.savefig(f'{out_prefix}_identity_percentage_histogram.pdf')
    plt.close()

    # 2) alignment length distribution (all)
    identity_min = merged_df_sorted['average_percent_identity'].min()
    identity_max = merged_df_sorted['average_percent_identity'].max()
    bins = np.linspace(identity_min, identity_max, num=21)
    merged_df_sorted['identity_bin'] = pd.cut(merged_df_sorted['average_percent_identity'], bins=bins, include_lowest=True)
    grouped_all = merged_df_sorted.groupby('identity_bin')['alignment_length_bp'].agg(['mean','std']).reset_index()
    bin_centers_all = [i.mid for i in grouped_all['identity_bin']]
    mean_lengths_all = grouped_all['mean']
    std_lengths_all = grouped_all['std']

    plt.figure(figsize=(10,6))
    plt.plot(bin_centers_all, mean_lengths_all, color='red', label='Mean Alignment Length')
    plt.fill_between(bin_centers_all,
                     mean_lengths_all - std_lengths_all,
                     mean_lengths_all + std_lengths_all,
                     color='red', alpha=0.3, label='Std Dev')
    plt.xlabel('Average Percent Identity')
    plt.ylabel('Alignment Length (bp)')
    plt.title('Alignment Length Distribution by Identity (All Data)')
    plt.legend()
    plt.savefig(f'{out_prefix}_alignment_length_distribution_all.pdf')
    plt.close()

    # 3) alignment length distribution (≥20kb)
    big_mask = merged_df_sorted['alignment_length_bp'] >= 20000
    merged_20kb = merged_df_sorted[big_mask].copy()
    merged_20kb['identity_bin'] = pd.cut(merged_20kb['average_percent_identity'], bins=bins, include_lowest=True)
    grouped_20kb = merged_20kb.groupby('identity_bin')['alignment_length_bp'].agg(['mean','std']).reset_index()
    bin_centers_20kb = [i.mid for i in grouped_20kb['identity_bin']]
    mean_lengths_20kb = grouped_20kb['mean']
    std_lengths_20kb = grouped_20kb['std']

    plt.figure(figsize=(10,6))
    plt.plot(bin_centers_20kb, mean_lengths_20kb, color='green', label='Mean Alignment Length (≥20kb)')
    plt.fill_between(bin_centers_20kb,
                     mean_lengths_20kb - std_lengths_20kb,
                     mean_lengths_20kb + std_lengths_20kb,
                     color='green', alpha=0.3, label='Std Dev')
    plt.xlabel('Average Percent Identity')
    plt.ylabel('Alignment Length (bp)')
    plt.title('Alignment Length Distribution by Identity (≥20kb)')
    plt.legend()
    plt.savefig(f'{out_prefix}_alignment_length_distribution_20kb.pdf')
    plt.close()

    # 4) Trend curve of alignment length by average identity
    num_bins = 30
    bins2 = np.linspace(identity_min, identity_max, num=num_bins+1)
    merged_df_sorted['identity_bin'] = pd.cut(merged_df_sorted['average_percent_identity'], bins=bins2, include_lowest=True)
    grouped_trend = merged_df_sorted.groupby('identity_bin')['alignment_length_bp'].agg(['mean','std']).reset_index()
    bin_centers = [i.mid for i in grouped_trend['identity_bin']]
    mean_alignment = grouped_trend['mean']
    std_alignment = grouped_trend['std']

    plt.figure(figsize=(10,6))
    plt.plot(bin_centers, mean_alignment, color='purple', label='Mean Alignment Length')
    plt.fill_between(bin_centers,
                     mean_alignment - std_alignment,
                     mean_alignment + std_alignment,
                     color='purple', alpha=0.3, label='Std Dev')
    plt.xlabel('Average Percent Identity (%)')
    plt.ylabel('Alignment Length (bp)')
    plt.title('Trend of Alignment Length by Average Percent Identity')
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.savefig(f'{out_prefix}_average_identity_vs_alignment_length_trend_curve.pdf')
    plt.close()

    # 5) KDE Plot with log Y-axis
    clean_df = merged_df_sorted.dropna(subset=['average_percent_identity','alignment_length_bp'])
    clean_df = clean_df[np.isfinite(clean_df['average_percent_identity']) & np.isfinite(clean_df['alignment_length_bp'])]
    clean_df = clean_df[clean_df['alignment_length_bp'] > 0]

    plt.figure(figsize=(10,6))
    kde = sns.kdeplot(
        x=clean_df['average_percent_identity'],
        y=clean_df['alignment_length_bp'],
        cmap='viridis',
        fill=True,
        thresh=0,
        levels=20
    )
    plt.yscale('log')
    plt.xlabel('Average Percent Identity (%)')
    plt.ylabel('Alignment Length (bp) [Log Scale]')
    plt.title('KDE of Average Percent Identity vs. Alignment Length (Log Scale)')
    plt.colorbar(kde.collections[0], label='Density')
    plt.savefig(f'{out_prefix}_average_identity_vs_alignment_length_kde_log.pdf')
    plt.close()

def main():
    parser = argparse.ArgumentParser(
        description="Single script that splits FASTA, creates BLAST DB, runs BLAST, merges synteny, produces plots."
    )
    parser.add_argument("--fasta", "-f", required=True, help="Path to the large FASTA file.")
    parser.add_argument("--output_dir", "-o", required=True, help="Output directory for chunked FASTA, BLAST DB, results, synteny, and plots.")
    parser.add_argument("--chunk_size", type=int, default=20000, help="Chunk size for splitting FASTA.")
    parser.add_argument("--threads", type=int, default=8, help="Number of threads for BLAST.")

    args = parser.parse_args()
    fasta_file = os.path.abspath(args.fasta)
    output_dir = os.path.abspath(args.output_dir)
    chunk_size = args.chunk_size
    threads = args.threads

    # Create main output dir if needed
    ensure_directory_exists(output_dir)

    # 1) Split FASTA
    chunked_fasta_dir = os.path.join(output_dir, "chunked_fasta")
    split_fasta(fasta_file, chunked_fasta_dir, chunk_size)

    # 2) Create BLAST DB
    db_path = os.path.join(output_dir, BLAST_DB_NAME)
    create_blast_db(chunked_fasta_dir, db_path)

    # 3) Run BLAST -> filtered lines -> concatenated
    temp_blast_dir = os.path.join(output_dir, TEMP_BLAST_DIRNAME)
    run_blast(chunks_dir=chunked_fasta_dir,
              db_path=db_path,
              temp_blast_dir=temp_blast_dir,
              threads=threads)

    # 4) Load final BLAST results
    concatenated_blast_path = os.path.join(temp_blast_dir, "concatenated_blast.txt")
    final_blast_txt = os.path.join(output_dir, BLAST_RESULTS_BASENAME)
    shutil.copy(concatenated_blast_path, final_blast_txt)  # Put a copy in main output
    processed_csv = os.path.join(output_dir, DETAILED_CSV_BASENAME)
    blast_df = load_and_process_blast_data(final_blast_txt, processed_csv)

    # 5) Merge synteny blocks
    synteny_df = process_synteny_blocks(blast_df)
    synteny_csv = os.path.join(output_dir, SYNTENY_CSV_BASENAME)
    synteny_df.to_csv(synteny_csv, index=False)

    # 6) Generate plots
    plot_prefix = os.path.join(output_dir, PLOT_PREFIX)
    generate_plots(synteny_df, plot_prefix)

    log("All steps completed successfully!")
    log(f"Final BLAST text: {final_blast_txt}")
    log(f"Detailed CSV: {processed_csv}")
    log(f"Synteny CSV: {synteny_csv}")
    log("Plots saved with prefix: " + plot_prefix)

if __name__ == "__main__":
    main()
