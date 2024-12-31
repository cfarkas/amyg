#!/usr/bin/env python3
import argparse
import os
import sys
import subprocess
import shutil
import logging
import time
import signal
import threading


###############################################################################
# GLOBAL SETTINGS
###############################################################################
REQUIRED_TOOLS = [
    "gffread",
    "gffcompare",
    "stringtie",
    "blastn",
    "gmap",
    "bedtools",
    "samtools",
    "TransDecoder.LongOrfs",
    "TransDecoder.Predict",
    "seqkit"
]

GREEN = "\033[92m"
RESET = "\033[0m"

logger = logging.getLogger("pipeline")
logger.setLevel(logging.DEBUG)

fh = logging.FileHandler("pipeline.log", mode="w")
fh.setLevel(logging.DEBUG)
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)

formatter = logging.Formatter(
    "%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S"
)
fh.setFormatter(formatter)
ch.setFormatter(formatter)

logger.addHandler(fh)
logger.addHandler(ch)


###############################################################################
# HELPER FUNCTIONS
###############################################################################
def log_green_info(message):
    """Print key steps in green to highlight them."""
    logger.info(f"{GREEN}{message}{RESET}")


def kill_docker_containers():
    """
    Kills any running containers based on myorg/annotate_env:latest image.
    Useful to ensure we don't leave stuck containers if user presses Ctrl+C.
    """
    logger.info("Killing all running containers from image myorg/annotate_env:latest (if any)...")
    cmd = "docker ps -q --filter=ancestor=myorg/annotate_env:latest | xargs -r docker kill || true"
    subprocess.run(cmd, shell=True)


def handle_sigint(signum, frame):
    """
    Handle Ctrl+C gracefully:
    - kill all Docker containers from myorg/annotate_env:latest
    - exit
    """
    logger.error("Ctrl+C caught! Terminating processes and exiting...")
    kill_docker_containers()
    sys.exit(1)


# Register the SIGINT handler
signal.signal(signal.SIGINT, handle_sigint)


def run_cmd(cmd, shell=True):
    """Run a shell command on the host and exit if it fails."""
    logger.debug(f"Running command: {cmd}")
    result = subprocess.run(cmd, shell=shell)
    if result.returncode != 0:
        logger.error(f"Command failed: {cmd}")
        sys.exit(1)


def run_pipeline_command(cmd, use_conda, use_docker, output_dir):
    """
    Run a pipeline command either in conda or docker environment.
    - If docker, mount output_dir as /data and run in container with `myorg/annotate_env:latest`.
    - If conda, run 'conda run -n annotate_env bash -c "cd {output_dir} && ..."'.
    - Otherwise, run on the host.
    """
    if use_docker:
        if not os.path.isabs(output_dir):
            logger.error("output_dir must be absolute when using docker.")
            sys.exit(1)

        uid = os.getuid()
        gid = os.getgid()
        full_cmd = (
            f"docker run --rm "
            f"-v {output_dir}:/data "
            f"-w /data "
            f"--user {uid}:{gid} "
            f"myorg/annotate_env:latest "
            f"bash -c \"{cmd}\""
        )
    elif use_conda:
        # In conda mode, we just run conda with local directory change
        full_cmd = f"conda run -n annotate_env bash -c \"cd {output_dir} && {cmd}\""
    else:
        full_cmd = cmd

    run_cmd(full_cmd)


###############################################################################
# BLAST progress monitor: track lines in transcriptome.swissprot
###############################################################################
class SwissprotMonitor(threading.Thread):
    """
    A background thread that, every 'interval' seconds,
    prints how many lines are in 'transcriptome.swissprot' if it exists.
    """
    def __init__(self, file_path, interval=60):
        super().__init__()
        self.file_path = file_path
        self.interval = interval
        self._stop_event = threading.Event()

    def run(self):
        while not self._stop_event.is_set():
            time.sleep(self.interval)
            if os.path.isfile(self.file_path):
                # count lines
                with open(self.file_path, 'r') as f:
                    n_lines = sum(1 for _ in f)
                logger.info(f"[BLAST progress] '{self.file_path}' has {n_lines} lines so far...")
            else:
                logger.info("[BLAST progress] transcriptome.swissprot not created yet...")

    def stop(self):
        self._stop_event.set()


def run_gawn_with_monitor(gawn_command, file_path, use_conda, use_docker, output_dir):
    """
    Runs GAWN in the foreground while a background thread
    prints the line count of 'transcriptome.swissprot' every 60s.
    """
    monitor = SwissprotMonitor(file_path=file_path, interval=60)
    monitor.start()

    try:
        # run pipeline command in foreground
        run_pipeline_command(gawn_command, use_conda, use_docker, output_dir)
    finally:
        # stop monitor
        monitor.stop()
        monitor.join()


###############################################################################
# ENVIRONMENT FILES (conda/docker)
###############################################################################
def environment_yml_content():
    """
    Return the environment.yml content, including:
    - procps-ng for 'ps',
    - parallel for parallel commands,
    - tqdm for progress bars in annotate_gtf.py.
    """
    return """\
name: annotate_env
channels:
  - bioconda
  - conda-forge
  - defaults
dependencies:
  - stringtie=2.2.1
  - gffcompare=0.11.2
  - gffread=0.12.7
  - blast=2.13.0
  - gmap=2021.08.25
  - bedtools=2.30.0
  - samtools=1.17
  - transdecoder=5.5.0
  - r-base=4.1.3
  - seqkit=2.3.1
  - parallel
  - procps-ng
  - tqdm
  - python=3.9
  - numpy
  - pandas
  - matplotlib
  - seaborn
  - biopython
  - intervaltree
"""


def dockerfile_content():
    """
    Return the Dockerfile content used to build the image with the conda environment.
    """
    return """\
FROM continuumio/miniconda3:4.8.2
WORKDIR /opt

COPY environment.yml /tmp/environment.yml
RUN conda env create -f /tmp/environment.yml && conda clean -afy

ENV PATH /opt/conda/envs/annotate_env/bin:$PATH
WORKDIR /data
"""


def write_environment_yml(content):
    with open("environment.yml", "w") as f:
        f.write(content)
    log_green_info("environment.yml written.")


def write_dockerfile(extra_files=None):
    lines = [
        "FROM continuumio/miniconda3:4.8.2",
        "WORKDIR /opt",
        "COPY environment.yml /tmp/environment.yml",
        "RUN conda env create -f /tmp/environment.yml && conda clean -afy",
        "ENV PATH /opt/conda/envs/annotate_env/bin:$PATH",
        "WORKDIR /data"
    ]
    if extra_files:
        for host_path, container_name in extra_files:
            lines.append(f"COPY {container_name} /data/{container_name}")

    with open("Dockerfile", "w") as f:
        f.write("\n".join(lines) + "\n")
    log_green_info("Dockerfile written (with extra files if any).")


###############################################################################
# Checking Tools
###############################################################################
def conda_run_which(tool):
    """Check if a tool is available in the 'annotate_env' conda environment."""
    cmd = f"conda run -n annotate_env which {tool}"
    result = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return (result.returncode == 0)


def docker_run_which(tool):
    """Check if a tool is available in the 'annotate_env' Docker image."""
    cmd = f"docker run --rm myorg/annotate_env:latest conda run -n annotate_env which {tool}"
    result = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return (result.returncode == 0)


def install_missing_tools_conda(missing):
    """Install any missing tools in the conda environment."""
    log_green_info("Installing missing tools in conda environment...")
    pkgs = " ".join(missing)
    cmd = f"conda run -n annotate_env conda install -c bioconda -c conda-forge -y {pkgs}"
    run_cmd(cmd)


def rebuild_docker_with_missing_tools(missing):
    """Add missing tools to environment.yml and rebuild the Docker image."""
    log_green_info("Updating environment.yml to include missing tools for Docker...")
    with open("environment.yml", "r") as f:
        lines = f.readlines()
    new_lines = []
    insert_index = None
    for i, line in enumerate(lines):
        new_lines.append(line)
        if line.strip().startswith("dependencies:"):
            insert_index = i + 1
    if insert_index is not None:
        for tool in missing:
            new_lines.insert(insert_index, f"  - {tool}\n")

    with open("environment.yml", "w") as f:
        f.writelines(new_lines)

    log_green_info("Rebuilding Docker image with updated environment.yml...")
    write_dockerfile()
    run_cmd("docker build . -t myorg/annotate_env:latest")


def verify_tools_conda():
    """Check each required tool in conda environment and install if missing."""
    log_green_info("Verifying tools in conda environment...")
    missing = []
    for tool in REQUIRED_TOOLS:
        logger.info(f"Checking tool: {tool}")
        if not conda_run_which(tool):
            logger.warning(f"{tool} not found in conda environment.")
            missing.append(tool)
    if missing:
        install_missing_tools_conda(missing)
        for tool in missing:
            if not conda_run_which(tool):
                logger.error(f"{tool} still not found after installation. Exiting.")
                sys.exit(1)
        log_green_info("All missing tools installed successfully in conda environment.")
    else:
        log_green_info("All required tools are present in the conda environment.")


def verify_tools_docker():
    """Check each required tool in Docker image and rebuild if missing."""
    log_green_info("Verifying tools in docker image...")
    missing = []
    for tool in REQUIRED_TOOLS:
        logger.info(f"Checking tool: {tool}")
        if not docker_run_which(tool):
            logger.warning(f"{tool} not found in Docker image.")
            missing.append(tool)
    if missing:
        log_green_info("Some tools are missing in Docker image. Attempting to rebuild image with missing tools.")
        rebuild_docker_with_missing_tools(missing)
        for tool in missing:
            if not docker_run_which(tool):
                logger.error(f"{tool} still not found after rebuild. Exiting.")
                sys.exit(1)
        log_green_info("All missing tools are now present in the Docker image.")
    else:
        log_green_info("All required tools are present in the Docker image.")


###############################################################################
# Env Installation
###############################################################################
def install_conda_env():
    """Create the conda environment and verify required tools."""
    log_green_info("Installing conda environment...")
    write_environment_yml(environment_yml_content())
    run_cmd("conda env create -f environment.yml")
    log_green_info("Conda environment 'annotate_env' created successfully.")
    verify_tools_conda()
    log_green_info("::: Installation Complete. Exiting. :::")
    sys.exit(0)


def install_docker_image():
    """Build the Docker image and verify required tools."""
    log_green_info("Installing docker image...")
    write_environment_yml(environment_yml_content())
    write_dockerfile()
    run_cmd("docker build . -t myorg/annotate_env:latest")
    log_green_info("Docker image 'myorg/annotate_env:latest' built successfully.")
    verify_tools_docker()
    log_green_info("::: Installation Complete. Exiting. :::")
    sys.exit(0)


###############################################################################
# Checking existence of conda env / docker image
###############################################################################
def conda_env_exists():
    """Check if the 'annotate_env' conda environment exists on the host."""
    cmd = "conda env list | grep annotate_env"
    result = subprocess.run(cmd, shell=True)
    return (result.returncode == 0)


def docker_image_exists():
    """Check if the 'myorg/annotate_env:latest' Docker image exists."""
    cmd = "docker images | grep myorg/annotate_env"
    result = subprocess.run(cmd, shell=True)
    return (result.returncode == 0)


###############################################################################
# Purge
###############################################################################
def purge_all_envs():
    """
    Remove the conda environment 'annotate_env' and the Docker image 'myorg/annotate_env:latest',
    then exit.
    """
    logger.info("Purging all environments (conda + docker)...")
    run_cmd("conda remove -n annotate_env --all -y || true")
    run_cmd("docker rmi myorg/annotate_env:latest -f || true")
    logger.info("All environments purged. Exiting now.")
    sys.exit(0)


###############################################################################
# Copy SwissProt for conda
###############################################################################
def copy_swissprot_conda(db_dir, gawn_dir):
    """
    In conda mode, copy SwissProt .pin, .psq, etc. from db_dir => gawn/03_data
    Then GAWN references them by local path e.g. '03_data/swissprot'
    """
    logger.info("Copying SwissProt DB into GAWN/03_data for conda mode to avoid /data references.")
    three_data = os.path.join(gawn_dir, "03_data")
    os.makedirs(three_data, exist_ok=True)
    run_cmd(f"cp -v {db_dir}/swissprot.* {three_data}/")


###############################################################################
# MAIN
###############################################################################
def main():
    log_green_info("Starting pipeline script...")

    parser = argparse.ArgumentParser(
        description="annotation pipeline that aims to annotate a de novo sequenced genomes (draft or complete) using RNA-seq evidence."
    )
    parser.add_argument("--install", choices=["conda", "docker"], help="Install environment and exit.")
    parser.add_argument("--use_conda", action="store_true", help="Run commands in conda env")
    parser.add_argument("--use_docker", action="store_true", help="Run commands in docker image")
    parser.add_argument("--threads", type=int, default=10, help="Number of CPUs (NCPUS) for gawn_config.sh")
    parser.add_argument("--force", action="store_true", help="Overwrite database and gawn_config.sh if present")
    parser.add_argument("--purge_all_envs", action="store_true", help="Remove the conda env and docker image, then exit.")
    parser.add_argument("-o", "--output", help="Output directory (must exist)")
    parser.add_argument("-a", help="StringTie GTF")
    parser.add_argument("-g", help="Reference genome (in fasta format)")

    args = parser.parse_args()

    # if purge
    if args.purge_all_envs:
        purge_all_envs()

    # if install
    if args.install == "conda":
        install_conda_env()
    elif args.install == "docker":
        install_docker_image()

    # normal mode
    if not args.install:
        if not args.a or not args.g or not args.output:
            logger.error("arguments -a, -g, --threads and -o must be provided when not using --install")
            sys.exit(1)

    use_conda = args.use_conda
    use_docker = args.use_docker
    threads = args.threads
    force = args.force
    a = args.a
    g = args.g
    output_dir = args.output

    # Convert to absolute, print
    if a:
        a = os.path.abspath(a)
        print(f"Full path for -a: {a}")
    if g:
        g = os.path.abspath(g)
        print(f"Full path for -g: {g}")
    if output_dir:
        output_dir = os.path.abspath(output_dir)
        print(f"Full path for -o: {output_dir}")

    if not args.install:
        # Check environment usage
        if use_conda and not conda_env_exists():
            logger.error("Conda environment 'annotate_env' not found. Run `--install conda` first.")
            sys.exit(1)
        if use_docker and not docker_image_exists():
            logger.error("Docker image 'myorg/annotate_env:latest' not found. Run `--install docker` first.")
            sys.exit(1)

        if not os.path.isdir(output_dir):
            logger.error(f"Output directory: {output_dir} not found. Please create it first.")
            sys.exit(1)

        db_dir = os.path.join(output_dir, "database")
        gawn_config_path = os.path.join(output_dir, "gawn_config.sh")

        # check if we need force
        if (os.path.exists(db_dir) or os.path.exists(gawn_config_path)) and not force:
            logger.error("Database or gawn_config.sh already exist in the output directory. Use --force to overwrite.")
            sys.exit(1)

        # force => remove
        if force:
            if os.path.exists(db_dir):
                shutil.rmtree(db_dir)
            if os.path.exists(gawn_config_path):
                os.remove(gawn_config_path)

        # (Optional) if args.minirun => slice GTF/FASTA (function omitted)
        # if args.minirun:
        #     a, g = create_mini_run(a, g, output_dir)

        # Download SwissProt
        log_green_info("Downloading and preparing SwissProt database inside output_dir/database...")
        os.makedirs(db_dir, exist_ok=True)
        run_cmd(f"wget -P {db_dir} ftp://ftp.ncbi.nlm.nih.gov/blast/db/swissprot.tar.gz")
        run_cmd(f"gunzip {os.path.join(db_dir,'swissprot.tar.gz')}")
        run_cmd(f"tar -xvf {os.path.join(db_dir,'swissprot.tar')} -C {db_dir}")

        # create initial gawn_config
        with open(gawn_config_path, "w") as gc:
            gc.write("#!/bin/bash\n")
            gc.write(f"NCPUS={threads}\n")
            gc.write("SKIP_GENOME_INDEXING=1\n")
            gc.write("GENOME_NAME=\"genome.fasta\"\n")
            gc.write("TRANSCRIPTOME_NAME=\"transcriptome.fasta\"\n")
            if use_docker:
                gc.write('SWISSPROT_DB="/data/database/swissprot"\n')
            else:
                gc.write('SWISSPROT_DB="TO_BE_REPLACED"\n')
            gc.write("#\n")

        # copy GTF/FASTA to output for environment
        a_filename = os.path.basename(a)
        g_filename = os.path.basename(g)
        shutil.copy(a, os.path.join(output_dir, a_filename))
        shutil.copy(g, os.path.join(output_dir, g_filename))

        # Step 1: gffread
        logger.info("::: Step 1: Obtaining Transcripts in FASTA format with gffread :::")
        run_pipeline_command(
            f"gffread -w transcripts.fa -g {g_filename} {a_filename}",
            use_conda, use_docker, output_dir
        )

        # Step 2: GAWN
        logger.info("::: Step 2: Setting up and running GAWN for gene annotation :::")
        GAWN_DIR = os.path.join(output_dir, "gawn")
        if os.path.isdir(GAWN_DIR):
            shutil.rmtree(GAWN_DIR)
        os.makedirs(GAWN_DIR, exist_ok=True)

        run_pipeline_command("git clone https://github.com/enormandeau/gawn.git gawn",
                             use_conda, use_docker, output_dir)

        # Copy transcripts/genome to GAWN data
        shutil.copy(os.path.join(output_dir, "transcripts.fa"), os.path.join(GAWN_DIR, "03_data", "transcriptome.fasta"))
        shutil.copy(os.path.join(output_dir, g_filename), os.path.join(GAWN_DIR, "03_data", "genome.fasta"))

        # If conda => copy SwissProt DB -> GAWN/03_data, fix config
        if use_conda:
            copy_swissprot_conda(db_dir, GAWN_DIR)
            # fix gawn_config: "SWISSPROT_DB=03_data/swissprot"
            with open(gawn_config_path, "r") as f:
                lines = f.readlines()
            new_lines = []
            for line in lines:
                if "SWISSPROT_DB" in line and "TO_BE_REPLACED" in line:
                    new_lines.append('SWISSPROT_DB="03_data/swissprot"\n')
                else:
                    new_lines.append(line)
            with open(gawn_config_path, "w") as f:
                f.writelines(new_lines)

        shutil.copy(gawn_config_path, os.path.join(GAWN_DIR, "02_infos", "gawn_config.sh"))

        # Step 2a: run GAWN with a BLAST line-count progress monitor
        logger.info("::: Running GAWN pipeline with BLAST line-count progress :::")
        swissprot_file = os.path.join(GAWN_DIR, "04_annotation", "transcriptome.swissprot")
        run_gawn_with_monitor("cd gawn && ./gawn 02_infos/gawn_config.sh",
                              swissprot_file,
                              use_conda,
                              use_docker,
                              output_dir)

        # Step 3: Extract transcriptome hits
        logger.info("::: Step 3: Extracting transcriptome hits :::")
        hits_path = os.path.join(GAWN_DIR, "04_annotation", "transcriptome.hits")
        swissprot_path = os.path.join(GAWN_DIR, "04_annotation", "transcriptome.swissprot")
        if not os.path.isfile(hits_path) or not os.path.isfile(swissprot_path):
            logger.error("transcriptome.hits or transcriptome.swissprot not found. GAWN may have failed.")
            sys.exit(9999)

        shutil.copy(swissprot_path, output_dir)
        shutil.copy(hits_path, output_dir)

        # Step 4: TransDecoder
        logger.info("::: Step 4: Predicting coding regions using TransDecoder :::")
        td_dir = os.path.join(output_dir, "transcripts.fa.transdecoder_dir")
        if os.path.exists(td_dir):
            shutil.rmtree(td_dir)

        run_pipeline_command("TransDecoder.LongOrfs -t transcripts.fa",
                             use_conda, use_docker, output_dir)
        run_pipeline_command("TransDecoder.Predict -t transcripts.fa",
                             use_conda, use_docker, output_dir)

        # Step 5: final_results
        logger.info("::: Step 5: Copying outputs to final_results :::")
        FINAL_RESULTS_DIR = os.path.join(output_dir, "final_results")
        os.makedirs(FINAL_RESULTS_DIR, exist_ok=True)
        files_to_move = [
            "transcripts.fa.transdecoder_dir/longest_orfs.gff3",
            "transcripts.fa.transdecoder_dir/longest_orfs.cds",
            "transcripts.fa.transdecoder_dir/longest_orfs.pep",
            "transcriptome.hits",
            "transcriptome.swissprot",
            "transcripts.fa"
        ]
        for f in files_to_move:
            src = os.path.join(output_dir, f)
            if os.path.isfile(src):
                shutil.move(src, FINAL_RESULTS_DIR)

        # Step 6: transdecoder_results
        logger.info("::: Step 6: Organizing remaining outputs into transdecoder_results :::")
        TRANSDECODER_RESULTS_DIR = os.path.join(output_dir, "transdecoder_results")
        os.makedirs(TRANSDECODER_RESULTS_DIR, exist_ok=True)

        for f in os.listdir(output_dir):
            if f.startswith("transcripts.fa.transdecoder."):
                shutil.move(os.path.join(output_dir, f), TRANSDECODER_RESULTS_DIR)

        for f in os.listdir(os.getcwd()):
            if f.startswith("pipeliner.") and f.endswith(".cmds"):
                shutil.move(os.path.join(os.getcwd(), f), TRANSDECODER_RESULTS_DIR)

        dirs_to_move = [
            "transcripts.fa.transdecoder_dir",
            "transcripts.fa.transdecoder_dir.__checkpoints",
            "transcripts.fa.transdecoder_dir.__checkpoints_longorfs"
        ]
        for d in dirs_to_move:
            dpath = os.path.join(os.getcwd(), d)
            if os.path.isdir(dpath):
                shutil.move(dpath, TRANSDECODER_RESULTS_DIR)

        logger.info("::: TransDecoder-related files moved to transdecoder_results directory :::")

        # Step 7: Annotate GTF
        logger.info("::: Step 7: Annotating GTF :::")
        logger.info("::: Downloading annotate_gtf.py script :::")
        run_pipeline_command(
            "curl -O https://raw.githubusercontent.com/cfarkas/annotate_my_genomes/master/additional_scripts/annotate_gtf.py",
            use_conda, use_docker, output_dir
        )
        run_pipeline_command("chmod +x annotate_gtf.py", use_conda, use_docker, output_dir)

        gtf_basename = os.path.basename(a)
        INPUT_GTF = gtf_basename
        HITS_FILE = "gawn/04_annotation/transcriptome.hits"
        ANNOTATION_TABLE = "gawn/05_results/transcriptome_annotation_table.tsv"
        TEMP_ANNOTATED_GTF = "final_annotated.gtf"

        # Check existence
        host_gtf_path = os.path.join(output_dir, gtf_basename)
        if not os.path.isfile(host_gtf_path):
            logger.error(f"Input GTF file not found in host output_dir: {host_gtf_path}")
            sys.exit(9999)

        host_hits_path = os.path.join(output_dir, HITS_FILE)
        if not os.path.isfile(host_hits_path):
            logger.error(f"Hits file not found in host output_dir: {host_hits_path}")
            sys.exit(9999)

        host_table_path = os.path.join(output_dir, ANNOTATION_TABLE)
        if not os.path.isfile(host_table_path):
            logger.error(f"Annotation table not found in host output_dir: {host_table_path}")
            sys.exit(9999)

        run_pipeline_command(
            f"python annotate_gtf.py {INPUT_GTF} {HITS_FILE} {ANNOTATION_TABLE} {TEMP_ANNOTATED_GTF}",
            use_conda,
            use_docker,
            output_dir
        )
        logger.info("::: GTF Annotation Completed :::")

        local_annotated_gtf = os.path.join(output_dir, "final_annotated.gtf")
        if os.path.isfile(local_annotated_gtf):
            shutil.move(local_annotated_gtf, FINAL_RESULTS_DIR)
            logger.info("::: Annotated GTF moved to final_results folder :::")
        else:
            logger.warning("No final_annotated.gtf found in host output_dir. Check if annotate_gtf.py ran correctly.")

        # Step 9: optional copying of transcriptome_annotation_table
        # table_src = os.path.join(output_dir, "gawn", "05_results", "transcriptome_annotation_table.tsv")
        # if os.path.isfile(table_src):
        #     shutil.copy(table_src, FINAL_RESULTS_DIR)
        #     logger.info("::: transcriptome_annotation_table.tsv copied to final_results :::")

        # Check if leftover
        contents = os.listdir(output_dir)
        exclude = {
            "final_results",
            "transdecoder_results",
            "database",
            "gawn_config.sh",
            "transcripts.fa",
            "gawn"
        }
        leftover = [c for c in contents if c not in exclude]
        if leftover:
            TIMESTAMP = time.strftime("%Y%m%d_%H%M%S")
            NEW_DIR = os.path.join(output_dir, f"annotate_my_genomes_{TIMESTAMP}")
            os.makedirs(NEW_DIR, exist_ok=True)
            logger.info(f"Moving organized results to {NEW_DIR} due to existing content in the output directory")
            if os.path.exists(FINAL_RESULTS_DIR):
                shutil.move(FINAL_RESULTS_DIR, os.path.join(NEW_DIR, "final_results"))
            if os.path.exists(os.path.join(output_dir, "transdecoder_results")):
                shutil.move(os.path.join(output_dir, "transdecoder_results"), os.path.join(NEW_DIR, "transdecoder_results"))
            FINAL_DIR = os.path.join(NEW_DIR, "final_results")
        else:
            FINAL_DIR = os.path.join(output_dir, "final_results")

        logger.info("::: Moving GAWN Directory :::")
        GAWN_DIR = os.path.join(output_dir, "gawn")
        if os.path.exists(GAWN_DIR):
            if 'NEW_DIR' in locals():
                logger.info(f"Moving GAWN directory to {NEW_DIR}")
                shutil.move(GAWN_DIR, NEW_DIR)
            else:
                logger.info(f"Moving GAWN directory to {FINAL_DIR}")
                shutil.move(GAWN_DIR, FINAL_DIR)

        logger.info("::: GAWN Directory Moved :::")

        logger.info("::: Navigating to Final Results Directory :::")
        # Step 8: Extracting CDS and peptide sequences
        logger.info("::: Step 8: Extracting CDS and peptide sequences to FASTA format :::")
        run_pipeline_command("gffread -g transcripts.fa -y cds.pep longest_orfs.gff3",
                             use_conda, use_docker, FINAL_DIR)
        run_pipeline_command("gffread -g transcripts.fa -x cds.fa longest_orfs.gff3",
                             use_conda, use_docker, FINAL_DIR)
        run_pipeline_command("gffread -g transcripts.fa -w three_prime_cds_five_prime.fa longest_orfs.gff3",
                             use_conda, use_docker, FINAL_DIR)

        logger.info("::: CDS and peptide sequences have been extracted to FASTA format :::")
        logger.info("::: Pipeline completed :::")


if __name__ == "__main__":
    main()