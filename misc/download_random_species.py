#!/usr/bin/env python3

import os
import random
import subprocess
import sys
import shutil
import argparse
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor, as_completed

ASSEMBLY_SUMMARY = "assembly_summary_refseq.txt"


def parse_args():
    parser = argparse.ArgumentParser(
        description="Randomly pick assemblies from a chosen group, ensure distinct species, "
                    "check SRA runs with specified library strategies, and download data. "
                    "If SRA downloads fail due to size or other errors, remove the species folder. "
                    "If genome + reads are already present, skip. You can also specify "
                    "a target species (e.g., 'Homo sapiens') to download instead of random picks."
    )
    parser.add_argument(
        "--group",
        choices=["all", "vertebrate", "invertebrate", "bacteria"],
        default="all",
        help="Which taxonomic group to filter on (default=all). Ignored if --target_species is set."
    )
    parser.add_argument(
        "--species", type=int, default=5,
        help="Number of distinct species to download (default=5). Ignored if --target_species is set."
    )
    parser.add_argument(
        "--target_species", type=str, default=None,
        help="If set, download only assemblies for this exact organism name (case-insensitive). "
             "Example: --target_species 'Homo sapiens'. This overrides the random approach."
    )
    parser.add_argument(
        "--seed", type=int, default=None,
        help="Random seed for reproducibility (default=None)."
    )
    parser.add_argument(
        "--min_sra", type=int, default=1,
        help="Minimum number of matching SRA runs required for a BioSample (default=1)."
    )
    parser.add_argument(
        "--max_sra", type=int, default=5,
        help="Maximum number of SRA runs to download per species (default=5)."
    )
    parser.add_argument(
        "--strategies", type=str, default="RNA-Seq",
        help="Comma-separated list of library strategies to require, e.g. 'RNA-Seq,Iso-Seq'. "
             "Default is 'RNA-Seq'."
    )
    parser.add_argument(
        "--threads", type=int, default=1,
        help="Number of parallel threads to use for SRA run downloads per species (default=1)."
    )
    return parser.parse_args()


def download_assembly_summary():
    """Download assembly_summary_refseq.txt if not present."""
    if not os.path.exists(ASSEMBLY_SUMMARY):
        url = ("ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/"
               "assembly_summary_refseq.txt")
        print(f"Downloading {ASSEMBLY_SUMMARY} ...")
        subprocess.run(["curl", "-O", url], check=True)
    else:
        print(f"{ASSEMBLY_SUMMARY} already exists, skipping download.")


def load_records_filtered_by_group(selected_group="all"):
    """
    Load all assemblies from 'assembly_summary_refseq.txt',
    filter by 'group' if not 'all'.

    group might be 'vertebrate_mammalian', 'vertebrate_other', 'invertebrate',
    'bacteria', 'plant', etc. We'll substring-match the user-chosen --group.

    Returns: list of dict with keys [gcf, biosample, taxid, species, group]
    """
    records = []
    with open(ASSEMBLY_SUMMARY, "r", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.strip().split("\t")
            if len(cols) < 25:
                continue

            # Columns (0-based):
            #  0 = assembly_accession (GCF_xxx)
            #  2 = biosample
            #  6 = species_taxid
            #  7 = organism_name
            # 24 = group
            gcf = cols[0]
            biosample = cols[2]
            taxid = cols[6]
            species_name = cols[7]
            group_value = cols[24].lower()  # e.g. 'vertebrate_mammalian', 'invertebrate'...

            if selected_group != "all":
                if selected_group not in group_value:
                    continue

            rec = {
                "gcf": gcf,
                "biosample": biosample,
                "taxid": taxid,
                "species": species_name,
                "group": group_value
            }
            records.append(rec)
    return records


def load_records_for_target_species(target_species):
    """
    Load all assemblies from 'assembly_summary_refseq.txt',
    and return only those whose organism_name matches target_species
    (case-insensitive exact match).

    Returns: list of dict with keys [gcf, biosample, taxid, species, group]
    """
    target_lower = target_species.strip().lower()
    matches = []
    with open(ASSEMBLY_SUMMARY, "r", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.strip().split("\t")
            if len(cols) < 25:
                continue
            gcf = cols[0]
            biosample = cols[2]
            taxid = cols[6]
            species_name = cols[7]
            group_value = cols[24].lower()

            if species_name.strip().lower() == target_lower:
                rec = {
                    "gcf": gcf,
                    "biosample": biosample,
                    "taxid": taxid,
                    "species": species_name,
                    "group": group_value
                }
                matches.append(rec)

    return matches


def build_sra_query(biosample, strategies):
    """
    Build an Entrez query to retrieve runs for the given BioSample
    that match any of the specified library strategies.
      e.g. SAMNxxxx AND (RNA-Seq[Strategy] OR Iso-Seq[Strategy])
    """
    if not biosample or biosample.lower() == "na":
        return None

    strategy_list = [s.strip() for s in strategies.split(",") if s.strip()]
    if not strategy_list:
        return biosample  # no strategy filter

    or_clauses = [f"{st}[Strategy]" for st in strategy_list]
    or_block = " OR ".join(or_clauses)
    return f"{biosample} AND ({or_block})"


def check_sra_data(biosample, strategies):
    """
    Use Entrez Direct (esearch+efetch) to find runs for a given BioSample.
    Return a list of run IDs (empty if none).
    If there's an error or 'PhraseNotFound', we treat it as zero results.
    """
    query = build_sra_query(biosample, strategies)
    if not query:
        return []

    cmd_esearch = ["esearch", "-db", "sra", "-query", query]
    cmd_efetch = ["efetch", "-format", "runinfo"]

    try:
        esearch_proc = subprocess.Popen(cmd_esearch, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)
        efetch_proc = subprocess.Popen(cmd_efetch, stdin=esearch_proc.stdout,
                                       stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)
        esearch_proc.stdout.close()
        output, _ = efetch_proc.communicate()
        esearch_proc.wait()
        efetch_proc.wait()
    except Exception:
        # If we fail for any reason, treat as no runs found
        return []

    lines = output.strip().split("\n")
    if len(lines) <= 1:
        return []  # no runs

    header = lines[0].split(",")
    if "Run" not in header:
        return []

    run_index = header.index("Run")
    runs = []
    for row in lines[1:]:
        parts = row.split(",")
        if len(parts) > run_index:
            run_acc = parts[run_index].strip()
            if run_acc:
                runs.append(run_acc)

    return runs


def download_assembly(gcf, species_folder):
    """Download genome (FASTA) + GTF + GFF3 via 'datasets' and unzip into species_folder."""
    zip_file = f"{gcf}.zip"
    cmd = [
        "datasets", "download", "genome", "accession", gcf,
        "--include", "genome",
        "--include", "gtf",
        "--include", "gff3",  # <--- ADDED LINE: also download GFF3
        "--filename", zip_file
    ]
    subprocess.run(cmd, check=True)

    # Unzip into species folder
    unzip_cmd = ["unzip", "-o", zip_file, "-d", species_folder]
    subprocess.run(unzip_cmd, check=True)


def _prefetch_run(run_acc, species_folder):
    """Helper function to run 'prefetch' on a single run."""
    # By default, prefetch has a 20GB limit. If you need bigger files:
    #   cmd = ["prefetch", "--max-size", "200GB", "--output-directory", species_folder, run_acc]
    cmd = ["prefetch", "--output-directory", species_folder, run_acc]
    subprocess.run(cmd, check=True)


def cleanup_species_folder(folder_path):
    """
    Remove the entire species folder if it exists.
    Called when downloads fail.
    """
    if os.path.isdir(folder_path):
        shutil.rmtree(folder_path, ignore_errors=True)


def already_has_genome_and_sra(species_folder):
    """
    Check if the folder has a genome/FASTA, at least one .sra file, AND a GFF3.
    If so, return True => skip.
    Criteria:
    - We look for *.fna, *.fa, *.fasta, or .gtf to confirm genome is present
    - We look for .sra files to confirm runs
    - We ALSO now look for .gff or .gff3 to confirm GFF3 is present
    """
    if not os.path.isdir(species_folder):
        return False

    has_sra = False
    has_genome = False
    has_gff = False

    for root, dirs, files in os.walk(species_folder):
        for f in files:
            f_lower = f.lower()
            if f_lower.endswith(".sra"):
                has_sra = True
            if (f_lower.endswith(".fna") or f_lower.endswith(".fa") or
                f_lower.endswith(".fasta") or f_lower.endswith(".gtf")):
                has_genome = True
            # GFF3 check
            if f_lower.endswith(".gff") or f_lower.endswith(".gff3"):
                has_gff = True

    # Now we only skip if we have all three: genome, SRA, GFF
    return (has_sra and has_genome and has_gff)


def download_sra_runs(run_ids, species_folder, max_sra=5, threads=1):
    """
    Download up to 'max_sra' runs with 'prefetch' into species_folder in parallel.
    """
    run_subset = run_ids[:max_sra]

    with ThreadPoolExecutor(max_workers=threads) as executor:
        futures = [executor.submit(_prefetch_run, r, species_folder) for r in run_subset]
        for _ in tqdm(as_completed(futures), total=len(futures), desc="Parallel SRA downloads", leave=False):
            pass


def main():
    args = parse_args()

    if args.seed is not None:
        random.seed(args.seed)

    # 1) Download assembly_summary if needed
    download_assembly_summary()

    # 2) If user specified --target_species, we load only that species
    if args.target_species:
        all_records = load_records_for_target_species(args.target_species)
        if not all_records:
            print(f"No assemblies found for species='{args.target_species}'. Exiting.")
            sys.exit(0)

        # For safety, let's treat it like we might get multiple assemblies, 
        # but we'll keep only distinct taxids so we don't re-download the same species multiple times.
        chosen_species = []
        chosen_taxids = set()

        for rec in all_records:
            if rec["taxid"] not in chosen_taxids:
                chosen_species.append(rec)
                chosen_taxids.add(rec["taxid"])

        print(f"Found {len(chosen_species)} assemblies for target_species='{args.target_species}'.")
    else:
        # Else, we do the random approach with filtering by group
        all_records = load_records_filtered_by_group(args.group)
        if not all_records:
            print(f"No assemblies found for group={args.group}. Exiting.")
            sys.exit(0)

        random.shuffle(all_records)

        chosen_species = []
        chosen_taxids = set()

        # TQDM bar for scanning
        print("Scanning assemblies to find species with sufficient runs ...")
        for rec in tqdm(all_records, desc="Scanning species", unit="assembly"):
            if len(chosen_species) >= args.species:
                break

            taxid = rec["taxid"]
            if taxid in chosen_taxids:
                continue

            run_ids = check_sra_data(rec["biosample"], args.strategies)
            if len(run_ids) < args.min_sra:
                continue

            chosen_species.append(rec)
            chosen_taxids.add(taxid)

        if len(chosen_species) < args.species:
            print(f"Warning: Could not find {args.species} distinct species with >= "
                  f"{args.min_sra} runs for strategies={args.strategies}. "
                  f"Found only {len(chosen_species)}.")

    # 3) Now download data for the chosen species
    for rec in tqdm(chosen_species, desc="Downloading data", unit="species"):
        gcf = rec["gcf"]
        species_name = rec["species"]
        biosample = rec["biosample"]

        species_folder = species_name.replace(" ", "_")
        run_ids = check_sra_data(biosample, args.strategies)

        # If we already have genome + SRA + GFF => skip
        if already_has_genome_and_sra(species_folder):
            print(f"Skipping {species_name} because genome + SRA + GFF3 are already present.")
            continue

        print(f"\nSpecies={species_name} (GCF={gcf}), Found {len(run_ids)} runs matching "
              f"strategies={args.strategies}. Downloading up to {args.max_sra} in parallel={args.threads}...")

        os.makedirs(species_folder, exist_ok=True)

        # Download assembly (FASTA + GTF + GFF3)
        try:
            download_assembly(gcf, species_folder)
        except subprocess.CalledProcessError as e:
            print(f"Error downloading assembly for {gcf}: {e}")
            # Cleanup folder
            cleanup_species_folder(species_folder)
            continue

        # Download SRA runs
        try:
            download_sra_runs(run_ids, species_folder, max_sra=args.max_sra, threads=args.threads)
        except subprocess.CalledProcessError as e:
            print(f"Error in parallel SRA downloads for species={species_name}: {e}")
            # Cleanup folder
            cleanup_species_folder(species_folder)
            continue

    print("Done.")


if __name__ == "__main__":
    # Check for required tools
    if shutil.which("datasets") is None:
        print("Error: 'datasets' CLI not found. Install from NCBI.")
        sys.exit(1)
    if shutil.which("esearch") is None or shutil.which("efetch") is None:
        print("Error: Entrez Direct (esearch, efetch) not found.")
        sys.exit(1)
    if shutil.which("prefetch") is None:
        print("Error: SRA Toolkit 'prefetch' not found.")
        sys.exit(1)

    main()
