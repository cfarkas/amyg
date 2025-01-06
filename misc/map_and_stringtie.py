#!/usr/bin/env python3

"""
map_and_stringtie.py (using GFF only)

1) For each run, produce a sorted BAM (no intermediate SAM) by piping hisat2 output to samtools sort,
   then samtools index it.
2) If --merge_runs is set, merge them into one final sorted BAM before StringTie.
3) fasterq-dump => (detect single vs. paired) => HISAT2 => samtools sort => samtools index => StringTie
4) pigz (parallel gzip) is used if available, else gzip, if --use_compression is set (forced with -f).

We assume the annotation is in a .gff (or .gff3) file, which we pass directly to StringTie.

Conda environment example:
    conda create -n amyg_benchmarking -c conda-forge -c bioconda \
        python=3.9 pigz sra-tools hisat2 samtools stringtie tqdm -y
    conda activate amyg_benchmarking

Usage:
    ./map_and_stringtie.py [options]

Required tools:
- pigz/gzip (if --use_compression)
- sra-tools (fasterq-dump)
- hisat2
- samtools
- stringtie
- python + tqdm
"""

import os
import sys
import argparse
import subprocess
import shutil
from glob import glob
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor, as_completed


def parse_args():
    parser = argparse.ArgumentParser(
        description="Recursively map reads with HISAT2 (single/paired autodetect) and call StringTie, "
                    "sorting & indexing all BAMs without writing intermediate SAM. "
                    "Now pigz/gzip uses force (-f) to avoid prompts."
    )
    parser.add_argument(
        "--threads", type=int, default=1,
        help="Number of CPU threads to use for alignment, stringtie, compression, etc. (default=1)."
    )
    parser.add_argument(
        "--working_dir", type=str, default=".",
        help="Base directory to search for species folders (default=current dir)."
    )
    parser.add_argument(
        "--merge_runs", action="store_true",
        help="If set, merge all BAMs from multiple runs into one final sorted BAM before StringTie."
    )
    parser.add_argument(
        "--verbose", action="store_true",
        help="Print extra debug info."
    )
    parser.add_argument(
        "--use_compression", action="store_true",
        help="If set, compress FASTQs after fasterq-dump using pigz if available, else gzip (forced with -f)."
    )
    # NEW ARGUMENT: --split_3
    parser.add_argument(
        "--split_3", action="store_true",
        help="Use 'fasterq-dump --split-3' instead of '--split-files' (helpful for partial paired data)."
    )
    return parser.parse_args()


def find_compressor():
    """Return 'pigz' if found, else 'gzip'."""
    if shutil.which("pigz"):
        return "pigz"
    else:
        return "gzip"


def run_cmd(cmd, cwd=None, verbose=False):
    """
    Helper to run a shell command with error checking.
    """
    if verbose:
        print(f"[CMD] {' '.join(cmd)} (cwd={cwd or '.'})")
    subprocess.run(cmd, check=True, cwd=cwd)


def find_species_folders(base_dir, verbose=False):
    """
    Recursively scan for directories containing:
      - ncbi_dataset/data/GCF_*/*.fna (genome)
      - SRR*/*.sra (reads)
    Return a list of such folder paths.
    """
    species_folders = []
    if verbose:
        print(f"[INFO] Scanning {base_dir} for species directories...")
    for root, dirs, files in os.walk(base_dir):
        fna = glob(os.path.join(root, "ncbi_dataset", "data", "GCF_*", "*.fna"))
        sra = glob(os.path.join(root, "SRR*", "*.sra"))
        if fna and sra:
            species_folders.append(root)
    species_folders = list(set(species_folders))
    if verbose:
        print(f"[INFO] Found {len(species_folders)} species folders.")
    return species_folders


def build_hisat2_index(genome_fna, index_prefix, threads=1, verbose=False):
    """
    Build HISAT2 index if not already present.
    """
    index_exists = all(os.path.exists(f"{index_prefix}.{i}.ht2") for i in range(1, 9))
    if index_exists:
        if verbose:
            print(f"[INFO] HISAT2 index at {index_prefix}.*.ht2 found, skipping build.")
        return
    cmd = ["hisat2-build", "-p", str(threads), genome_fna, index_prefix]
    if verbose:
        print(f"[INFO] Building HISAT2 index for {genome_fna}")
    run_cmd(cmd, verbose=verbose)


def run_fasterq_dump(sra_path, threads=1, use_compression=False, verbose=False, split_3=False):
    """
    Convert .sra -> FASTQ using either:
      - fasterq-dump --split-files  (default), or
      - fasterq-dump --split-3      (if split_3=True)

    Return a list of FASTQ(s). If single-end, 1 file; if paired, 2 files.
    If use_compression, also compress them with pigz/gzip (forced with -f).
    """
    sra_dir = os.path.dirname(sra_path)
    sra_base = os.path.splitext(os.path.basename(sra_path))[0]
    if verbose:
        print(f"[INFO] Running fasterq-dump on {sra_path}")

    # 1) Build the command
    if split_3:
        # user asked for --split-3
        cmd = ["fasterq-dump", "--split-3", "--threads", str(threads), sra_path]
    else:
        cmd = ["fasterq-dump", "--split-files", "--threads", str(threads), sra_path]

    run_cmd(cmd, cwd=sra_dir, verbose=verbose)

    # 2) Identify produced FASTQs
    produced_fastqs = []
    single = os.path.join(sra_dir, f"{sra_base}.fastq")
    r1 = os.path.join(sra_dir, f"{sra_base}_1.fastq")
    r2 = os.path.join(sra_dir, f"{sra_base}_2.fastq")

    # If --split-3, partial pairs might appear in the single file <run>.fastq
    # We keep the same single vs. paired detection.
    if os.path.exists(single):
        produced_fastqs.append(single)

    if os.path.exists(r1):
        produced_fastqs.append(r1)
    if os.path.exists(r2):
        produced_fastqs.append(r2)

    # 3) Compress if requested (force with -f)
    if use_compression and produced_fastqs:
        comp = find_compressor()
        if verbose:
            print(f"[INFO] Using {comp} to compress {len(produced_fastqs)} FASTQ(s) with -f.")
        for fq in produced_fastqs:
            if comp == "pigz":
                run_cmd([comp, "-p", str(threads), "-f", fq], verbose=verbose)
            else:
                run_cmd([comp, "-f", fq], verbose=verbose)
        produced_fastqs = [fq + ".gz" for fq in produced_fastqs]

    return produced_fastqs


def run_hisat2_alignment(index_prefix, fastqs, out_prefix, threads=1, verbose=False):
    """
    Directly pipe HISAT2 -> samtools sort -> final .bam.
    Return path to final sorted + indexed .bam.

    Detect if single-end (fastqs=1) or paired-end (fastqs=2).
    """

    if not fastqs:
        if verbose:
            print("[WARN] No FASTQs to align. Skipping alignment.")
        return None

    sorted_bam = out_prefix + ".sorted.bam"

    import subprocess
    if len(fastqs) == 1:
        # single-end
        fq = fastqs[0]
        if verbose:
            print(f"[INFO] Aligning single-end: {fq} => {sorted_bam} (pipe).")
        hisat2_cmd = [
            "hisat2",
            "-p", str(threads),
            "-x", index_prefix,
            "-U", fq,
            "--quiet"
        ]
    elif len(fastqs) == 2:
        # paired-end
        fq1, fq2 = fastqs
        if verbose:
            print(f"[INFO] Aligning paired-end: {fq1}, {fq2} => {sorted_bam} (pipe).")
        hisat2_cmd = [
            "hisat2",
            "-p", str(threads),
            "-x", index_prefix,
            "-1", fq1,
            "-2", fq2,
            "--quiet"
        ]
    else:
        if verbose:
            print(f"[WARN] Found {len(fastqs)} FASTQs (expected 1 or 2). Skipping alignment.")
        return None

    if verbose:
        print("[INFO] Launching hisat2, piping directly to samtools sort.")

    hisat2_proc = subprocess.Popen(
        hisat2_cmd,
        stdout=subprocess.PIPE,
        stderr=None,
        text=False
    )
    sort_cmd = ["samtools", "sort", "-@", str(threads), "-o", sorted_bam, "-"]
    sort_proc = subprocess.Popen(
        sort_cmd,
        stdin=hisat2_proc.stdout,
        stderr=None
    )
    hisat2_proc.stdout.close()

    ret1 = hisat2_proc.wait()
    ret2 = sort_proc.wait()
    if ret1 != 0:
        raise subprocess.CalledProcessError(ret1, hisat2_cmd)
    if ret2 != 0:
        raise subprocess.CalledProcessError(ret2, sort_cmd)

    # index
    run_cmd(["samtools", "index", sorted_bam], verbose=verbose)

    return sorted_bam


def merge_and_sort_bams(bam_list, out_prefix, threads=1, verbose=False):
    """
    Merge multiple BAMs => produce merged.unsorted.bam => sort => merged.sorted.bam => index => remove unsorted
    Return final merged.sorted.bam
    """
    merged_unsorted = out_prefix + ".unsorted.bam"
    merged_sorted = out_prefix + ".sorted.bam"

    merge_cmd = ["samtools", "merge", "-@", str(threads), merged_unsorted] + bam_list
    run_cmd(merge_cmd, verbose=verbose)

    sort_cmd = ["samtools", "sort", "-@", str(threads), "-o", merged_sorted, merged_unsorted]
    run_cmd(sort_cmd, verbose=verbose)

    index_cmd = ["samtools", "index", merged_sorted]
    run_cmd(index_cmd, verbose=verbose)

    if os.path.exists(merged_unsorted):
        os.remove(merged_unsorted)

    return merged_sorted


def run_stringtie(in_bam, annotation_gff, out_gtf, threads=1, verbose=False):
    """
    stringtie -p <threads> -L -v -a 4 -G annotation_gff -o out_gtf in_bam
    """
    cmd = [
        "stringtie",
        "-p", str(threads),
        "-L",
        "-v",
        "-a", "4",
        "-G", annotation_gff,
        "-o", out_gtf,
        in_bam
    ]
    if verbose:
        print(f"[INFO] Running StringTie on {in_bam} using GFF={annotation_gff}")
    run_cmd(cmd, verbose=verbose)


def main():
    args = parse_args()
    base_dir = os.path.abspath(args.working_dir)
    threads = args.threads
    verbose = args.verbose

    if verbose:
        print(f"[INFO] Searching for species in {base_dir} ...")
    species_folders = find_species_folders(base_dir, verbose=verbose)
    if not species_folders:
        print("[ERROR] No species folders found. Exiting.")
        sys.exit(1)

    for species_dir in tqdm(species_folders, desc="Species", unit="species"):
        if verbose:
            print(f"\n[INFO] Processing species: {species_dir}")

        # find genome + GFF
        fna_list = glob(os.path.join(species_dir, "ncbi_dataset", "data", "GCF_*", "*.fna"))
        gff_list = glob(os.path.join(species_dir, "ncbi_dataset", "data", "GCF_*", "*.gff*"))
        if not fna_list or not gff_list:
            if verbose:
                print(f"[WARN] Missing .fna or .gff in {species_dir}, skipping.")
            continue

        genome_fna = fna_list[0]
        annotation_gff = gff_list[0]

        # Build HISAT2 index if needed
        index_prefix = os.path.join(species_dir, "hisat2_index")
        build_hisat2_index(genome_fna, index_prefix, threads=threads, verbose=verbose)

        # For each SRR folder
        sra_folders = glob(os.path.join(species_dir, "SRR*"))
        if not sra_folders:
            if verbose:
                print(f"[WARN] No SRR subfolders in {species_dir}.")
            continue

        all_bams = []
        for sra_fold in tqdm(sra_folders, desc="Runs", leave=False):
            sra_files = glob(os.path.join(sra_fold, "*.sra"))
            if not sra_files:
                if verbose:
                    print(f"[WARN] {sra_fold} has no .sra, skipping.")
                continue

            sra_path = sra_files[0]
            if verbose:
                print(f"[INFO] Converting SRA -> FASTQ -> alignment -> sorting -> indexing: {sra_path}")

            # Convert to FASTQ using --split-files or --split-3
            fastqs = run_fasterq_dump(
                sra_path,
                threads=threads,
                use_compression=args.use_compression,
                verbose=verbose,
                split_3=args.split_3  # <--- pass user choice
            )

            # Align => sorted BAM (via piping)
            out_prefix = os.path.join(sra_fold, "aligned")
            sorted_bam = run_hisat2_alignment(index_prefix, fastqs, out_prefix, threads=threads, verbose=verbose)
            if sorted_bam and os.path.exists(sorted_bam):
                all_bams.append(sorted_bam)

        if not all_bams:
            if verbose:
                print(f"[WARN] No BAM files for {species_dir}, skipping StringTie.")
            continue

        # Merge if requested
        if len(all_bams) == 1 or not args.merge_runs:
            final_bam = all_bams[0]
        else:
            merged_prefix = os.path.join(species_dir, "merged_runs")
            final_bam = merge_and_sort_bams(all_bams, merged_prefix, threads=threads, verbose=verbose)

        # Run StringTie
        out_gtf = os.path.join(species_dir, "transcripts.gtf")
        run_stringtie(final_bam, annotation_gff, out_gtf, threads=threads, verbose=verbose)
        if verbose:
            print(f"[INFO] Created {out_gtf}")

    print("[INFO] Pipeline complete!")


if __name__ == "__main__":
    # Check required tools
    required_tools = ["fasterq-dump", "hisat2", "samtools", "stringtie"]
    for tool in required_tools:
        if shutil.which(tool) is None:
            print(f"Error: '{tool}' not found in PATH. Install or activate environment, then retry.")
            sys.exit(1)

    main()
