# **amyg**: A Pipeline for De Novo Genomic Annotation of Non-Model Organisms

**amyg.py** is a Python-based annotation pipeline that aims to annotate a de novo sequenced genomes (draft or complete) using RNA-seq evidence. Currently the pipeline:
- Performs GTF processing from StringTie outputs  
- Generates gene annotation using [GAWN](https://github.com/enormandeau/gawn) with SwissProt/BLAST integration  
- Resolve transcriptome coding potential with **TransDecoder**, producing **longest ORFs**, **CDS**, and **peptide** sequences for each transcript.     

Currently, the pipeline can run through:
1. **Conda**  (an environment called `annotate_env` will be created in your system)
2. **Docker** (with an auto-built image `myorg/annotate_env:latest`)

## Synopsis
```
python3 amyg.py --help

usage: amyg.py [-h] [--install {conda,docker}] [--use_conda] [--use_docker] [--threads THREADS] [--force] [--purge_all_envs] [--dups]
               [--chunk_size CHUNK_SIZE] [-o OUTPUT] [-a A] [-g G]

annotation pipeline that aims to annotate a de novo sequenced genome using RNA-seq plus optional synteny BLAST for duplicates.

options:
  -h, --help            show this help message and exit
  --install {conda,docker}
                        Install environment and exit.
  --use_conda           Run commands in conda env
  --use_docker          Run commands in docker image
  --threads THREADS     Number of CPUs (NCPUs) for gawn_config.sh
  --force               Overwrite database and gawn_config.sh if present
  --purge_all_envs      Remove the conda env and docker image, then exit.
  --dups                Enable chunk-based synteny BLAST to find duplicates (will run amyg_syntenyblast.py).
  --chunk_size CHUNK_SIZE
                        Chunk size for synteny-based duplication step (only used if --dups is enabled).
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
- Docker image takes ~47.8 min to build in Ubuntu 24.04.1 LTS. We aimed to create a reproducible and robust local Docker image. Apologies for the delay. 

---

## Run
Currently there are two ways to run the pipeline:

### 1) Docker Mode
```
mkdir test_docker
python3 amyg.py \
  -a /path/to/my_genome.gtf \
  -g /path/to/my_genome.fasta \
  -o ./test_docker \
  --threads 25 \
  --use_docker \
  --force
```

- ```--threads 25``` sets number of cpus (NCPUs) for BLAST-based GAWN annotation.
- The output is placed in ```./test_docker```.

### 2) Conda Mode
```
mkdir test_conda
python3 amyg.py \
  -a /path/to/my_genome.gtf \
  -g /path/to/my_genome.fasta \
  -o ./test_conda \
  --threads 25 \
  --use_conda \
  --force
```
- The output is placed in ```./test_conda```.

#### Notes:

- **Ctrl+C** kills all running Docker containers, ensuring no stuck processes.
- ```--force``` overwrites existing database/ and gawn_config.sh if they are in the output folder. We reccomend to run the pipeline fresh using this flag. 

---

## Detailed Steps

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
--dups enables chunk-based synteny BLAST via amyg_syntenyblast.py to identify potential duplicated regions.
--chunk_size controls the size of each FASTA split for BLAST runs when --dups is used.

**Organizes** final results in `final_results/` subfolder and leftover TransDecoder outputs in `transdecoder_results/`.

## Interested in genome-wide duplications? please run the ```--dups``` flag

### 1) Docker Mode
```
mkdir test_docker
python3 amyg.py \
  -a /path/to/my_genome.gtf \
  -g /path/to/my_genome.fasta \
  -o ./test_docker \
  --threads 25 \
  --use_docker \
  --force \
  --dups
```
### 2) Conda Mode
```
mkdir test_conda
python3 amyg.py \
  -a /path/to/my_genome.gtf \
  -g /path/to/my_genome.fasta \
  -o ./test_conda \
  --threads 25 \
  --use_conda \
  --force \
  --dups
```
- Enabling ```--dups``` flag will also enable ```--chunk_size``` that will slice the genome (default at 20000 bp) and will test synteny comparing all fragments vs all, and at the end will reconstruct genomic segment with strong duplication evidence across the genome. Also, it will produce ```final_annotated_dups.gtf```which contains the annotation of duplicated genes on the ```final_annotated.gtf``` file    
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
