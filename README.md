![My Image](amyg.png)
# **amyg**: A Pipeline for De Novo Genomic Annotation of Non-Model Organisms compatible with Single Cell RNA-seq

**amyg** is a Python-based annotation pipeline that aims to annotate a de novo sequenced genomes (draft or complete) using RNA-seq evidence. Currently the pipeline:
- Performs GTF processing from StringTie outputs  
- Generates gene annotation using [GAWN](https://github.com/enormandeau/gawn) with SwissProt/BLAST integration  
- Resolve transcriptome coding potential with **TransDecoder**, producing **longest ORFs**, **CDS**, and **peptide** sequences for each transcript.     
- In single-cell mode, reconstruct GTF files from cell-type-specific bams (output of SComatic ```SplitBamCellTypes.py```), recognize novel genes and annotate.

Currently, the pipeline can run through:
1. **Conda**  (an environment called `annotate_env` will be created in your system)
2. **Docker** (with an auto-built image `myorg/annotate_env:latest`)

- See https://pypi.org/project/amyg/0.1.6/ 

## Synopsis
```
amyg --help

usage: amyg [-h] [--install {conda,docker}] [--use_conda] [--use_docker] [--threads THREADS] [--force] [--purge_all_envs] [--dups]
            [--chunk_size CHUNK_SIZE] [-o OUTPUT] [-a A] [-g G] [--egap_gff EGAP_GFF] [--single_cell] [--input_dir INPUT_DIR] [--ref_gtf REF_GTF]
            [--ref_fa REF_FA] [--overlap_frac OVERLAP_FRAC] [--preprocessing] [--bam BAM] [--continue]

annotation pipeline with optional single_cell or coverage-based preprocessing.

options:
  -h, --help            show this help message and exit
  --install {conda,docker}
                        Install environment and exit.
  --use_conda           Use conda env.
  --use_docker          Use docker image.
  --threads THREADS
  --force
  --purge_all_envs
  --dups
  --chunk_size CHUNK_SIZE
  -o OUTPUT, --output OUTPUT
                        Output directory
  -a A                  GTF from StringTie or param tuning
  -g G                  Reference genome (FASTA)
  --egap_gff EGAP_GFF   EGAP GFF for merging or reference checks
  --single_cell         If set => multi-bam single-cell logic => exit after done.
  --input_dir INPUT_DIR
                        Directory with .bam for single_cell
  --ref_gtf REF_GTF     Known ref GTF for single_cell mode
  --ref_fa REF_FA       Reference genome for single_cell final pass
  --overlap_frac OVERLAP_FRAC
                        Overlap fraction for single_cell mode.
  --preprocessing       If set => coverage-based param sets => pick best => override -a.
  --bam BAM             BAM used for coverage detection in preprocessing
  --continue            If set, continue the pipeline after preprocessing instead of exiting.
```

**amyg** is the next version of [annotate_my_genomes](https://github.com/cfarkas/annotate_my_genomes) but streamlines the installation and there is no need for separate config files. 

---

## Installation

Via pip: 
```
pip install amyg
```
Then, users can decide to install all requirements via conda or docker as follows: 

```bash
# 1) Install conda environment:
amyg --install conda

# 2) Install docker image:
amyg --install docker

# 3) Uninstall and purge old envs (optional):
amyg --purge_all_envs
```
- While Conda is faster, Docker image takes ~47.8 min to build in Ubuntu 24.04.1 LTS. We aimed to create a reproducible and robust local Docker image. Apologies for the delay. 

---

## Run
Currently there are two ways to run the pipeline:

### 1) Docker Mode
```
mkdir test_docker
amyg \
  -a /path/to/my_genome.gtf \
  -g /path/to/my_genome.fasta \
  -o ./test_docker \
  --threads 25 \
  --use_docker \
  --force
```

- ```--threads 25``` sets number of cpus (NCPUs) for BLAST-based GAWN annotation.
- The output is placed in ```./test_docker```. The main results of the pipeline will be inside i.e: ```./test_docker/amyg_20250101_150629/final_results/```

### 2) Conda Mode
```
mkdir test_conda
amyg \
  -a /path/to/my_genome.gtf \
  -g /path/to/my_genome.fasta \
  -o ./test_conda \
  --threads 25 \
  --use_conda \
  --force
```
- The output is placed in ```./test_conda```. The main results of the pipeline will be inside i.e: ```./test_conda/amyg_20250101_150629/final_results/```

#### Notes:

- **Ctrl+C** kills all running Docker containers, ensuring no stuck processes.
- ```--force``` overwrites existing database/ and gawn_config.sh if they are in the output folder. We reccomend to run the pipeline fresh using this flag. 

---

## Detailed Steps

1. **Download SwissProt**  
   - Automatically fetches `swissprot.tar.gz` from the NCBI BLAST FTP server and unpacks it into the `database/` folder.

2. **Create `gawn_config.sh`**  
   - **Docker mode** sets `SWISSPROT_DB` to `/data/database/swissprot`.  
   - **Conda mode** copies SwissProt into `gawn/03_data` and sets `SWISSPROT_DB` to `03_data/swissprot`.

3. **Run GAWN**  
   - BLAST progress is monitored every 60 seconds, logging how many lines appear in `transcriptome.swissprot`.

4. **TransDecoder**  
   - Discovers **longest ORFs** and **predicts coding regions**.

5. **Annotate GTF**  
   - Downloads `annotate_gtf.py` and merges final hits into `final_annotated.gtf`.
   - Outputs organized to `final_results/`, with any remaining TransDecoder files moved to `transdecoder_results/`.

6. **Usage with Optional `--dups`**
--dups enables chunk-based synteny BLAST via ```amyg_syntenyblast.py``` to identify potential duplicated regions.
--chunk_size controls the size of each FASTA split for BLAST runs when --dups is used.

**Organizes** final results in `final_results/` subfolder and leftover TransDecoder outputs in `transdecoder_results/`.

---
## Preprocessing your stringtie.gtf using ```--preprocessing``` flag

Sometimes you need to clean or prepare your GTF file before running the main annotation pipeline. The ```--preprocessing``` flag lets you do just that. Here's what it does in detail:

1) **(Recommended) Run ```unique_gene_id.py```**

This script ensures all gene_id fields in your GTF are truly unique. For any conflicting gene_ids (e.g., multiple transcripts with the same gene_name), it automatically appends a suffix to avoid collisions.
The output is a new GTF (e.g., ```mygtf.gtf``` => ```mygtf.unique_gene_id.gtf```).

2) **(Optional, but recommended) Run ```merge_stringtie_names.py``` with ```--egap_gff``` (NCBI Eukaryotic Genome Annotation Pipeline gff)**

If you have the NCBI Eukaryotic Genome Annotation Pipeline gff of your genome and provide --egap_gff ```/path/to/genome.gff```, the pipeline automatically also downloads and runs the ```merge_stringtie_names.py``` script.
That script further refines your GTF by comparing it to the reference (the “EGAP GFF”), ensuring consistent naming of transcripts and unifying gene_id vs. gene_name across transcripts and exons.
The final result is a new file (by default named ```annotated_and_renamed.gtf```), which can then be used in the main amyg pipeline.

#### 1) If you just want to unique‐ify your gene IDs:
```
amyg --preprocessing -a /path/to/mygtf.gtf
```

#### 2) If you want to unique‐ify and also want to merge your GTF with an EGAP reference for consistent naming
```
amyg --preprocessing \
  -a /path/to/mygtf.gtf \
  --egap_gff /path/to/genomic.gff
```
```transcripts_named.gtf``` will be produced, that can be input for amyg pipeline.

---
## Interested in genome-wide duplications? please run with ```--dups``` flag

### 1) Docker Mode
```
mkdir test_docker
amyg \
  -a /path/to/my_genome.gtf \
  -g /path/to/my_genome.fasta \
  -o ./test_docker \
  --threads 25 \
  --use_docker \
  --force \
  --dups \
  --chunk_size 20000
```
### 2) Conda Mode
```
mkdir test_conda
amyg \
  -a /path/to/my_genome.gtf \
  -g /path/to/my_genome.fasta \
  -o ./test_conda \
  --threads 25 \
  --use_conda \
  --force \
  --dups \
  --chunk_size 20000
```
- Enabling ```--dups``` flag will also enable ```--chunk_size``` that will slice the genome (default at 20000 bp) and will test synteny comparing all fragments vs all, and at the end will reconstruct genomic segment with strong duplication evidence across the genome. Also, it will produce ```final_annotated_dups.gtf```which contains the annotation of duplicated genes on the ```final_annotated.gtf``` file    
- The results of the pipeline will be inside i.e: ```./output_folder/amyg_20250101_150629/final_results/```
---

## Plot duplications
- Inside i.e.: ```amyg_20250101_150629/final_results``` users can do the following 
```
wget https://raw.githubusercontent.com/cfarkas/amyg/refs/heads/main/scripts/plot_dups.py
chmod 755 plot_dups.py

python plot_dups.py \
    -a transcriptome_annotation_table.tsv \
    -g final_annotated_dups.gtf \
    -s synteny_blocks.csv \
    -o ./
```
This script will produce two stacked bar plots. 

Plot A:
  - Sorted contigs with percentage of Ancient (green), Recent (orange), Other (gray) genes

Plot B:
  - Classify *all duplications* (either “ancient” or “recent”) as:
      * "intra-only" (blue)   if contig is found ONLY in self-synteny blocks
      * "inter-only" (red)    if contig is found ONLY in cross-synteny blocks
      * "both"       (black)  if contig is found in both self- & cross-synteny
      * "other"      (gray)   if duplication_type ∉ {ancient,recent} 
                              OR contig not found in synteny at all
    So the bar shows how we partition the entire set of duplicated genes 
    among (intra-only, inter-only, both). 

---
---
## Applying amyg to single-cell data
In addition to annotate de novo sequenced genomes, amyg is able to recognize and annotate novel genes from single-cell data. The input files for running amyg in single-cell mode are cell-type-specific bams (output from SComatic ```SplitBamCellTypes.py```), with which amyg will automatically reconstruct the GTF files (using stringtie). Then, it runs the main pipeline (download SwissProt, create gawn_config.sh, run GAWN, run TransDecoder, and annotate GTF). Consider that to run amyg in this mode, cell annotation must already have been performed.

To run amyg with single-cell mode, it is necessary to add the flag ```--single_cell``` and the reference fasta and gtf file (```--ref_fa``` and ```--ref_gtf```, respectively).

### 1) Docker Mode
```
mkdir test_docker_scMode
amyg \
  --input_dir /path/to/SplitBamCellTypes.bam \
  -g /path/to/my_genome.fasta \
  --single_cell \
  --ref_fa /path/to/my_genome.fasta \
  --ref_gtf /path/to/my_genome.gtf \
  -o ./test_docker_scMode\
  --threads 25 \
  --use_docker \
  --force
```
### 2) Conda Mode
```
mkdir test_conda_scMode
amyg \
  --input_dir /path/to/SplitBamCellTypes.bam \
  -g /path/to/my_genome.fasta \
  --single_cell \
  --ref_fa /path/to/my_genome.fasta \
  --ref_gtf /path/to/my_genome.gtf \
  -o ./test_conda_scMode \
  --threads 25 \
  --use_conda \
  --force 
```

Similar to the option without single-cell mode, the main results of the pipeline will be in the ```final_results/``` subfolder. The output files are:
- final_annotated-gtf
- longest_orfs.cds
- longest_orfs.gff3
- longest_orfs.pep
- transcriptome_annotation_table.tsv
- transcriptome.hits
- transcriptome.swissprot
- transcripts.fa

---
## Applying amyg to single-cell data
Additionally, to annotate de novo sequenced genomes, amyg is able to recognize and annotate novel genes from single-cell data. The input files for running amyg in single-cell mode are cell-type-specific bams (output from SComatic SplitBamCellTypes.py), with which amyg will automatically reconstruct the GTF files (using stringtie). Then, it runs the main pipeline (download SwissProt, create gawn_config.sh, run GAWN, run TransDecoder, and annotate GTF). 
To run amyg with single-cell mode is necessary to add the flag ```--single_cell``` and the reference fasta and gtf file (```--ref_fa``` and ```--ref_gtf``` respectively).

### 1) Docker Mode
```
mkdir test_docker_scMode
amyg \
  --input_dir /path/to/SplitBamCellTypes.bam \
  -g /path/to/my_genome.fasta \
  --single_cell \
  --ref_fa /path/to/my_genome.fasta \
  --ref_gtf /path/to/my_genome.gtf \
  -o ./test_docker_scMode\
  --threads 25 \
  --use_docker \
  --force
```
### 2) Conda Mode
```
mkdir test_conda_scMode
amyg \
  --input_dir /path/to/SplitBamCellTypes.bam \
  -g /path/to/my_genome.fasta \
  --single_cell \
  --ref_fa /path/to/my_genome.fasta \
  --ref_gtf /path/to/my_genome.gtf \
  -o ./test_conda_scMode \
  --threads 25 \
  --use_conda \
  --force 
```

Similar to the option without single-cell mode, the main results of the pipeline will be in the ```final_results/``` subfolder. The output files are:
- final_annotated-gtf
- longest_orfs.cds
- longest_orfs.gff3
- longest_orfs.pep
- transcriptome_annotation_table.tsv
- transcriptome.hits
- transcriptome.swissprot
- transcripts.fa

---
  
## Requirements

- **Python 3.7+**  
- **Miniconda** or **Docker** installed on your system  
- Enough disk space for BLAST DB and GTF/FASTA inputs

---

### Troubleshooting

**Ctrl+C** in the middle of a run 
Kills Docker containers so you don’t have to manually do it.

**Permission**  
Make sure you have write access to your output directory and local Docker permissions.

---

### License

This project is licensed under the MIT License.
