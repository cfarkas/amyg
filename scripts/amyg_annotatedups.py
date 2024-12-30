import sys
import pandas as pd
from pybedtools import BedTool
from tqdm import tqdm
import re
import os
import logging

# -------------------------------
# Step 1: Define Input and Output Files
# -------------------------------

# Check if the correct number of arguments is provided
if len(sys.argv) != 5:
    print("Usage: python amyg_annotatedups.py <gtf_file> <synteny_blocks.csv> <output_gtf_file> <log_file>")
    sys.exit(1)

gtf_file = sys.argv[1]               # GTF file
synteny_blocks_file = sys.argv[2]    # Synteny blocks CSV file
output_gtf = sys.argv[3]             # Output annotated GTF file
log_file = sys.argv[4]               # Log file for debugging

# -------------------------------
# Step 2: Setup Logging
# -------------------------------

logging.basicConfig(
    filename=log_file,
    filemode='w',
    level=logging.DEBUG,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

# -------------------------------
# Step 3: Load and Inspect Synteny Blocks
# -------------------------------

logging.info("Loading synteny_blocks.csv")
try:
    synteny_df = pd.read_csv(synteny_blocks_file)
    logging.info(f"synteny_blocks.csv loaded with {len(synteny_df)} rows.")
except Exception as e:
    logging.error(f"Error loading synteny_blocks.csv: {e}")
    sys.exit(1)

# Verify required columns
required_columns = [
    'query_contig', 'query_seq_start', 'query_seq_end',
    'subject_contig', 'subject_seq_start', 'subject_seq_end',
    'duplication_type'
]
missing_columns = set(required_columns) - set(synteny_df.columns)
if missing_columns:
    logging.error(f"Missing columns in synteny_blocks.csv: {missing_columns}")
    sys.exit(1)
else:
    logging.info("All required columns are present in synteny_blocks.csv.")

# -------------------------------
# Step 4: Verify Chromosome/Contig Naming Consistency
# -------------------------------

logging.info("Extracting unique contigs from synteny_blocks.csv")
unique_synteny_contigs = set(synteny_df['query_contig'].unique()).union(set(synteny_df['subject_contig'].unique()))
logging.info(f"Unique contigs in synteny_blocks.csv: {unique_synteny_contigs}")

def extract_unique_chromosomes(gtf_path):
    """Extracts unique chromosome/contig names from the GTF file."""
    contigs = set()
    with open(gtf_path, 'r') as gtf:
        for line in gtf:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 1:
                continue
            contigs.add(parts[0])
    return contigs

logging.info("Extracting unique contigs from genes.gtf")
unique_gtf_contigs = extract_unique_chromosomes(gtf_file)
logging.info(f"Unique contigs in genes.gtf: {unique_gtf_contigs}")

# Check for discrepancies
missing_in_gtf = unique_synteny_contigs - unique_gtf_contigs
missing_in_synteny = unique_gtf_contigs - unique_synteny_contigs

if missing_in_gtf:
    logging.warning(f"Contigs present in synteny_blocks.csv but missing in genes.gtf: {missing_in_gtf}")
if missing_in_synteny:
    logging.warning(f"Contigs present in genes.gtf but missing in synteny_blocks.csv: {missing_in_synteny}")

if missing_in_gtf or missing_in_synteny:
    logging.warning("Potential chromosome/contig name mismatches detected. Proceeding with caution.")

# -------------------------------
# Step 5: Convert Duplication Events to BedTool
# -------------------------------

logging.info("Preparing duplication regions for BedTool conversion.")

# Prepare query duplications
query_duplications = synteny_df[['query_contig', 'query_seq_start', 'query_seq_end', 'duplication_type']].copy()
query_duplications.columns = ['chrom', 'start', 'end', 'dup_type']

# Prepare subject duplications
subject_duplications = synteny_df[['subject_contig', 'subject_seq_start', 'subject_seq_end', 'duplication_type']].copy()
subject_duplications.columns = ['chrom', 'start', 'end', 'dup_type']

# Combine query and subject duplications
all_duplications = pd.concat([query_duplications, subject_duplications], ignore_index=True)

# Remove duplicates
before_drop = len(all_duplications)
all_duplications.drop_duplicates(inplace=True)
after_drop = len(all_duplications)
logging.info(f"Removed {before_drop - after_drop} duplicate duplication regions.")

# Check if duplication regions are present
if all_duplications.empty:
    logging.error("No duplication regions found after processing synteny_blocks.csv.")
    sys.exit(1)
else:
    logging.info(f"Total duplication regions to process: {len(all_duplications)}")

# Convert to BedTool
try:
    dup_bedtool = BedTool.from_dataframe(all_duplications)
    logging.info("Duplication regions converted to BedTool object successfully.")
except Exception as e:
    logging.error(f"Error converting duplication regions to BedTool: {e}")
    sys.exit(1)

# -------------------------------
# Step 6: Convert GTF Transcripts to BedTool
# -------------------------------

logging.info("Extracting transcript features from genes.gtf")

def extract_transcripts_from_gtf(gtf_path):
    """Extract transcript lines from GTF as a Pandas DataFrame (for BedTool)."""
    transcripts = []
    with open(gtf_path, 'r') as gtf:
        for line in gtf:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            if parts[2] != 'transcript':
                continue
            chrom = parts[0]
            try:
                start = int(parts[3]) - 1  # GTF is 1-based; BED is 0-based
                end = int(parts[4])       # BED end is exclusive
            except ValueError:
                logging.warning(f"Invalid coordinates in line: {line.strip()}")
                continue
            strand = parts[6]
            attr = parts[8]
            # Extract transcript_id and gene_id
            transcript_id_match = re.search(r'transcript_id "([^"]+)"', attr)
            gene_id_match = re.search(r'gene_id "([^"]+)"', attr)
            transcript_id = transcript_id_match.group(1) if transcript_id_match else 'NA'
            gene_id = gene_id_match.group(1) if gene_id_match else 'NA'
            transcripts.append([chrom, start, end, transcript_id, 0, strand, gene_id])
    df = pd.DataFrame(transcripts, columns=['chrom', 'start', 'end', 'transcript_id', 'score', 'strand', 'gene_id'])
    return df

transcripts_df = extract_transcripts_from_gtf(gtf_file)
logging.info(f"Extracted {len(transcripts_df)} transcript features from genes.gtf.")

if transcripts_df.empty:
    logging.error("No transcript features found in genes.gtf.")
    sys.exit(1)

try:
    transcripts_bedtool = BedTool.from_dataframe(transcripts_df)
    logging.info("Transcript features converted to BedTool object successfully.")
except Exception as e:
    logging.error(f"Error converting transcript features to BedTool: {e}")
    sys.exit(1)

# -------------------------------
# Step 7: Perform Overlap Analysis
# -------------------------------

logging.info("Performing overlap between transcripts and duplication regions.")

try:
    overlaps = transcripts_bedtool.intersect(dup_bedtool, wa=True, wb=True)
    overlap_data = []
    for entry in overlaps:
        # transcript fields
        transcript_id = entry[3]
        transcript_gene_id = entry[6]
        # duplication fields
        dup_type = entry[10]
        
        overlap_data.append({
            'transcript_id': transcript_id,
            'gene_id': transcript_gene_id,
            'duplication_type': dup_type
        })
    
    overlap_df = pd.DataFrame(overlap_data)
    logging.info(f"Found {len(overlap_df)} overlapping transcript-duplication entries.")
except Exception as e:
    logging.error(f"Error during overlap analysis: {e}")
    sys.exit(1)

# -------------------------------
# Step 8: Assign Duplication Types to Transcripts
# -------------------------------

if not overlap_df.empty:
    duplication_priority = {'recent': 2, 'ancient': 1}
    
    def determine_dup_type(dup_types):
        """Given a list of duplication types, pick the highest priority one."""
        sorted_types = sorted(dup_types, key=lambda d: duplication_priority.get(d, 0), reverse=True)
        return sorted_types[0]
    
    duplication_annotations = overlap_df.groupby('transcript_id')['duplication_type'].apply(list).reset_index()
    duplication_annotations['assigned_dup_type'] = duplication_annotations['duplication_type'].apply(determine_dup_type)
    logging.info("Assigned duplication types based on overlap analysis.")
else:
    duplication_annotations = pd.DataFrame(columns=['transcript_id', 'duplication_type', 'assigned_dup_type'])
    logging.warning("No overlaps found. All transcripts will be annotated as 'no_duplication'.")

# Merge duplication annotations with transcripts_df
transcripts_annotated = transcripts_df.merge(
    duplication_annotations[['transcript_id', 'assigned_dup_type']], on='transcript_id', how='left'
)
transcripts_annotated['assigned_dup_type'].fillna('no_duplication', inplace=True)

dup_type_counts = transcripts_annotated['assigned_dup_type'].value_counts()
logging.info(f"Duplication Type Distribution:\n{dup_type_counts}")

# -------------------------------
# Step 9: Annotate GTF with Duplication Types
# -------------------------------

logging.info("Annotating genes.gtf with duplication_type.")

dup_dict = pd.Series(
    transcripts_annotated.assigned_dup_type.values,
    index=transcripts_annotated.transcript_id
).to_dict()

def annotate_gtf_with_duplication(gtf_path, dup_mapping, output_path):
    """Annotate transcript lines of the GTF with duplication_type."""
    with open(gtf_path, 'r') as infile, open(output_path, 'w') as outfile:
        for line in tqdm(infile, desc="Annotating GTF"):
            if line.startswith('#'):
                outfile.write(line)
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                outfile.write(line)
                continue
            if parts[2] != 'transcript':
                outfile.write(line)
                continue
            attr = parts[8]
            match = re.search(r'transcript_id "([^"]+)"', attr)
            if match:
                transcript_id = match.group(1)
                dup_type = dup_mapping.get(transcript_id, 'no_duplication')
                if not attr.endswith(';'):
                    attr += ';'
                attr += f' duplication_type "{dup_type}";'
                parts[8] = attr
                modified_line = '\t'.join(parts) + '\n'
                outfile.write(modified_line)
            else:
                outfile.write(line)

try:
    annotate_gtf_with_duplication(gtf_file, dup_dict, output_gtf)
    logging.info(f"Annotation completed successfully. Output saved to {output_gtf}")
except Exception as e:
    logging.error(f"Error during GTF annotation: {e}")
    sys.exit(1)

# -------------------------------
# Step 10: Final Checks
# -------------------------------

try:
    with open(output_gtf, 'r') as annotated_gtf:
        duplication_found = False
        for line in annotated_gtf:
            if line.startswith('#'):
                continue
            if 'duplication_type' in line:
                duplication_found = True
                break
        if duplication_found:
            logging.info("duplication_type annotation verified in the output GTF file.")
        else:
            logging.warning("duplication_type annotation not found in the output GTF file.")
except Exception as e:
    logging.error(f"Error during final verification: {e}")
    sys.exit(1)
