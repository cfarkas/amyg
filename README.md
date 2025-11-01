![My Image](amyg.png)
# **amyg**: A Pipeline for De Novo Genomic Annotation of Non‑Model Organisms — with splice‑aware exon correction and single‑cell RNA‑seq support

**amyg** is a Python-based annotation pipeline to annotate de novo genomes (draft or complete) using RNA‑seq evidence. This version adds an optional **exon‑correction stage** (Portcullis + PASA or Mikado), **keeps single‑cell support**, and **removes** the previous StringTie **parameter‑tuning** loop and the **duplication plotting** steps.

**Core capabilities**
- Normalize and prepare **GTF/GFF3** (unique IDs; optional EGAP naming).
- Optional **splice‑aware exon correction** using **Portcullis** junctions, followed by **PASA** *or* **Mikado**.
- Generate functional annotation with **[GAWN](https://github.com/enormandeau/gawn)** (SwissProt/BLAST/DIAMOND integration).
- Resolve coding potential with **TransDecoder** (longest ORFs, CDS, peptide sequences).
- **Single‑cell mode**: build per‑BAM assemblies with StringTie, merge/reference‑name, harvest truly novel transcripts, and continue into the main pipeline.

You can run **amyg** via:
1. **Conda** (creates an environment named `annotate_env`), or
2. **Docker** (builds an image `myorg/annotate_env:latest`).

See also PyPI: https://pypi.org/project/amyg/

---

## Synopsis

```text
amyg --help

usage: amyg [-h] [--install {conda,docker}] [--use_conda] [--use_docker]
            [--threads THREADS] [--force] [--purge_all_envs]
            [-o OUTPUT] [-a A] [-g G] [--egap_gff EGAP_GFF]
            [--single_cell] [--input_dir INPUT_DIR] [--ref_gtf REF_GTF]
            [--ref_fa REF_FA] [--overlap_frac OVERLAP_FRAC]
            [--preprocessing]
            [--exon_correction {none,pasa,mikado}] [--bam BAM]

annotation pipeline with optional single_cell and exon-correction (Portcullis + PASA/Mikado).

options:
  -h, --help            show this help message and exit
  --install {conda,docker}
                        Install environment and exit.
  --use_conda           Use conda env.
  --use_docker          Use docker image.
  --threads THREADS     Number of threads.
  --force               Overwrite existing SwissProt DB / gawn_config in output.
  --purge_all_envs      Remove the conda env and docker image, then exit.
  -o OUTPUT, --output OUTPUT
                        Output directory (absolute path required with --use_docker).
  -a A                  Input annotation (GTF or GFF3).
  -g G                  Reference genome (FASTA).
  --egap_gff EGAP_GFF   EGAP GFF for naming/consistency during preprocessing.
  --single_cell         Enable single-cell mode (per-BAM assemblies + merge).
  --input_dir INPUT_DIR Directory containing per-cell/cluster BAMs.
  --ref_gtf REF_GTF     Reference GTF for single-cell naming/merge.
  --ref_fa REF_FA       Reference FASTA for single-cell final pass.
  --overlap_frac OVERLAP_FRAC
                        Overlap fraction for single-cell filtering (default: 0.05).
  --preprocessing       Normalize IDs and optional EGAP naming (no StringTie tuning).
  --exon_correction {none,pasa,mikado}
                        Splice-aware exon correction (requires --bam and -g if not 'none').
  --bam BAM             RNA-seq BAM for Portcullis (coordinate-sorted and indexed).
```

**Note:** Preprocessing now runs inline and **the pipeline continues automatically**. There is no `--continue` flag in this version.

---

## Installation

Install the package:
```bash
pip install amyg
```

Then choose an execution backend:

```bash
# 1) Create the conda environment:
amyg --install conda

# 2) Build the docker image:
amyg --install docker

# 3) (Optional) Purge both envs:
amyg --purge_all_envs
```

> **Conda vs Docker:** Conda is often quicker to set up; Docker offers a fully reproducible runtime. With Docker, pass an **absolute** `-o/--output` path.

---

## Typical runs

### A) Minimal run (no exon correction)
**Conda**
```bash
amyg \
  -a assembly.gtf \
  -g genome.fa \
  -o ./out_conda \
  --threads 16 \
  --use_conda
```

**Docker** (absolute output path)
```bash
amyg \
  -a /data/assembly.gtf \
  -g /data/genome.fa \
  -o /ABSOLUTE/path/out_docker \
  --threads 16 \
  --use_docker
```

Outputs land in `out_*/final_results/` (see _Outputs_ below).

---

### B) Preprocessing GTF with EGAP naming (no exon correction)
```bash
amyg \
  --preprocessing \
  -a assembly.gtf \
  -g genome.fa \
  --egap_gff egap.gff3 \
  -o ./out_preproc \
  --use_conda
```
This normalizes IDs and merges StringTie naming against EGAP; the run then continues into GAWN + TransDecoder automatically.

---

### C) Exon correction: **Portcullis + Mikado** (recommended when you have RNA‑seq BAM)
```bash
amyg \
  --preprocessing \
  --exon_correction mikado \
  --bam rnaseq.sorted.bam \
  -a assembly.gtf \
  -g genome.fa \
  --egap_gff egap.gff3 \
  --threads 16 \
  -o ./out_mikado \
  --use_conda
```
What happens: Portcullis computes high‑confidence junctions; Mikado selects corrected models, which then feed into GAWN + TransDecoder.

---

### D) Exon correction: **Portcullis + PASA (SQLite)**
```bash
amyg \
  --preprocessing \
  --exon_correction pasa \
  --bam rnaseq.sorted.bam \
  -a assembly.gff3 \
  -g genome.fa \
  --egap_gff egap.gff3 \
  --threads 16 \
  -o ./out_pasa \
  --use_conda
```
What happens: transcripts are extracted and aligned with GMAP via PASA; PASA‑refined models feed into GAWN + TransDecoder.

> **Tip:** BAM must be coordinate‑sorted and indexed (`.bai`). Genome FASTA should be indexed (`samtools faidx`).

---

### E) Single‑cell mode
Provide per‑cell/cluster BAMs plus reference GTF/FA; `amyg` will run StringTie per BAM, merge, name, filter **truly novel** transcripts, **then continue** into the main pipeline.

```bash
amyg \
  --single_cell \
  --input_dir ./cells_bam \
  --ref_gtf ref_annotation.gtf \
  --ref_fa  genome.fa \
  --threads 16 \
  -o ./out_sc \
  --use_conda
```

The single‑cell stage writes to `out_sc/amyg_singlecell_YYYYMMDD/` and flows into GAWN + TransDecoder automatically.

---

## Detailed steps

1. **(Optional) Preprocessing**  
   - Normalize IDs (via `unique_gene_id.py`).  
   - If `--egap_gff` is provided, merge StringTie naming against EGAP for consistent IDs.

2. **(Optional) Exon correction**  
   - **Portcullis** computes reliable splice junctions from `--bam`.  
   - Choose either:  
     - **Mikado**: configure → prepare → serialise (with junctions/ORFs) → pick best loci; or  
     - **PASA (SQLite)**: transcript alignment & model refinement using GMAP.  
   - Corrected models are converted to GTF and replace `-a` for downstream steps.

3. **SwissProt/GAWN**  
   - Downloads/unpacks SwissProt DB under `output/database/`, prepares `gawn_config.sh`, and runs **GAWN** (BLAST/DIAMOND).

4. **TransDecoder**  
   - Finds longest ORFs and predicts coding regions.

5. **Annotate GTF**  
   - Integrates SwissProt hits/annotation table into **`final_annotated.gtf`**.

6. **Packaging**  
   - Assembles deliverables into `final_results/` and archives other artifacts in `transdecoder_results/` and method‑specific subfolders.

---

## Outputs

Inside `OUTPUT/final_results/` you will typically find:
- `final_annotated.gtf`
- `transcripts.fa`
- `transcriptome.swissprot`
- `transcriptome.hits`
- `transcriptome_annotation_table.tsv`
- `longest_orfs.gff3`
- `longest_orfs.cds`
- `longest_orfs.pep`

Method‑specific folders (if used):
- `portcullis_out/` — junction evidence
- `mikado/` — corrected loci (`mikado.loci.gff3` + converted GTF)
- PASA SQLite files/config/logs in the output root

---

## Requirements

- **Python 3.9+**
- **Miniconda** or **Docker**
- Enough disk space for SwissProt/BLAST databases and intermediate files

---

## Troubleshooting

- **Docker output path** must be absolute (`-o /ABS/path/out`).  
- **BAM indexing**: `samtools sort -o rnaseq.sorted.bam raw.bam && samtools index rnaseq.sorted.bam`.  
- **Interrupting runs**: Ctrl+C cleanly stops Docker containers created by the script.  
- **Permissions**: ensure you can write to the output folder and have permission to run Docker (if used).

---

## License

MIT License.
