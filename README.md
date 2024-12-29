# **amyg**: A Pipeline for De Novo Genomic Annotation of Non-Model Organisms

**amyg.py** is a Python-based annotation pipeline that aims to annotate a de novo sequenced genomes (draft or complete) using RNA-seq evidence. **amyg.py**:
- Performs GTF processing from StringTie outputs  
- Generates gene annotation using [GAWN](https://github.com/enormandeau/gawn) with SwissProt/BLAST integration  
- Resolve transcriptome coding potential with **TransDecoder**, producing **longest ORFs**, **CDS**, and **peptide** sequences for each transcript.     

Currently, it supports two main environments:

1. **Conda**  
2. **Docker** (with an auto-built image `myorg/annotate_env:latest`)

## Synopsis
```
python3 amyg.py --help

usage: amyg.py [-h] [--install {conda,docker}] [--use_conda] [--use_docker] [--threads THREADS] [--force] [--purge_all_envs] [-o OUTPUT] [-a A] [-g G]

Run pipeline with environment, database, dynamic gawn_config.sh. Overwrite with --force.

options:
  -h, --help            show this help message and exit
  --install {conda,docker}
                        Install environment and exit.
  --use_conda           Run commands in conda env
  --use_docker          Run commands in docker image
  --threads THREADS     Number of CPUs (NCPUS) for gawn_config.sh
  --force               Overwrite database and gawn_config.sh if present
  --purge_all_envs      Remove the conda env and docker image, then exit.
  -o OUTPUT, --output OUTPUT
                        Output directory (must exist)
  -a A                  StringTie GTF
  -g G                  Reference genome (in fasta format)
```

**amyg** is the next version of [annotate_my_genomes](https://github.com/cfarkas/annotate_my_genomes) but streamlines the installation and there is no need for separate config files.

---

## Installation

```bash
# 1) Install conda environment:
python3 amyg.py --install conda

# 2) Install docker image:
python3 amyg.py --install docker

# 3) Uninstall and purge old envs (optional):
python3 amyg.py --purge_all_envs
```

---

## Run
Currently there are two ways to run the pipeline:

### 1) Docker Mode
```
mkdir test_docker
python3 amyg.py \
  -a /path/to/genes.gtf \
  -g /path/to/genome.fasta \
  --threads 25 \
  -o /absolute/path/to/test_docker \
  --use_docker
```
- ```--threads 25``` sets NCPUs in GAWN.
- The output is placed in /absolute/path/to/test_docker.
- **Ctrl+C** kills all running Docker containers, ensuring no stuck processes.

### 2) Conda Mode
```
mkdir test_conda
python3 amyg.py \
  -a gtf.unique_gene_id.gtf \
  -g genome.fasta \
  --threads 25 \
  -o ./test_conda \
  --use_conda \
  --force
```
- ```--force``` overwrites existing database/ and gawn_config.sh if they are in the output folder.
- The pipeline runs entirely within your local annotate_env conda environment.

---

### Detailed Steps

**1. Download SwissProt**  
The script downloads `swissprot.tar.gz` from the NCBI BLAST ftp server, unpacks it into `database/`.

**2. Create `gawn_config.sh`**  
- In Docker mode, `SWISSPROT_DB` is set to `/data/database/swissprot`.
- In Conda mode, the pipeline copies SwissProt locally into `gawn/03_data` and sets `SWISSPROT_DB` to `03_data/swissprot`.

**3. Runs GAWN**  
BLAST progress is monitored every 60s, printing how many lines `transcriptome.swissprot` has so far.

**3. TransDecoder**  
Identifies **longest ORFs** → Predicts **coding regions**.

**4. Annotates GTF**  
Downloads `annotate_gtf.py` and merges final hits into `final_annotated.gtf`.

**Organizes** final results in `final_results/` subfolder and leftover TransDecoder outputs in `transdecoder_results/`.

---

### Requirements

- **Python 3.7+**  
- **Miniconda** or **Docker** installed on your system  
- Enough disk space for BLAST DB and GTF/FASTA inputs

---

### Troubleshooting

**Ctrl+C**  
Kills Docker containers so you don’t have to manually do it.

**Permission**  
Make sure you have write access to your output directory and local Docker permissions.

---

### License

This project is licensed under the MIT License.



