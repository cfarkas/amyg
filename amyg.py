#!/usr/bin/env python3
# amyg.py
# - Keeps preprocessing (-a GTF) and --egap_gff merge/naming path.
# - Adds Portcullis + (PASA | Mikado) exon-correction options.
# - Keeps single-cell mode.
# - Updates conda/docker recipes to include portcullis, pasa, mikado.
#
# Notes:
#   * Exon-correction is optional and controlled with:  --exon_correction {none,pasa,mikado}
#   * Portcullis is always run when exon_correction != none (needs --bam and -g).
#   * PASA path writes a SQLite config and runs Launch_PASA_pipeline.pl (GMAP-based).
#   * Mikado path runs: configure → prepare → (optional) serialise → pick using Portcullis junctions.
#   * Outputs from PASA/Mikado are converted to GTF and replace -a for downstream steps.

import argparse
import os
import sys
import subprocess
import shutil
import logging
import time
import signal
import threading
import glob
import statistics
import datetime

###############################################################################
# GLOBAL SETTINGS
###############################################################################
REQUIRED_TOOLS = [
    # core IO / comparison
    "gffread",
    "gffcompare",
    "bedtools",
    "samtools",
    "seqkit",

    # transcript assembly / evidence
    "stringtie",                 # used in single-cell mode
    "blastn",
    "gmap",

    # coding sequences
    "TransDecoder.LongOrfs",
    "TransDecoder.Predict",

    # exon-correction options
    "portcullis",
    "mikado",
    "Launch_PASA_pipeline.pl",
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
    Kills any running containers from myorg/annotate_env:latest.
    """
    logger.info("Killing all running containers from image myorg/annotate_env:latest (if any)...")
    cmd = "docker ps -q --filter=ancestor=myorg/annotate_env:latest | xargs -r docker kill || true"
    subprocess.run(cmd, shell=True)


def handle_sigint(signum, frame):
    """
    Handle Ctrl+C: kill Docker containers, exit.
    """
    logger.error("Ctrl+C caught! Terminating processes and exiting...")
    kill_docker_containers()
    sys.exit(1)

signal.signal(signal.SIGINT, handle_sigint)


def run_cmd(cmd, shell=True):
    """Run a shell command on the host. Exit if it fails."""
    logger.debug(f"Running command: {cmd}")
    result = subprocess.run(cmd, shell=shell)
    if result.returncode != 0:
        logger.error(f"Command failed: {cmd}")
        sys.exit(1)


def run_pipeline_command(cmd, use_conda, use_docker, output_dir):
    """
    Run 'cmd' in either conda or docker environment. If docker, mount 'output_dir'.
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
            f"bash -lc \"{cmd}\""
        )
    elif use_conda:
        full_cmd = f"conda run -n annotate_env bash -lc \"cd {output_dir} && {cmd}\""
    else:
        full_cmd = cmd
    run_cmd(full_cmd)

###############################################################################
# SwissprotMonitor + run_gawn_with_monitor
###############################################################################
class SwissprotMonitor(threading.Thread):
    """
    Prints how many lines in transcriptome.swissprot every 'interval' seconds.
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
                with open(self.file_path, 'r') as f:
                    n_lines = sum(1 for _ in f)
                logger.info(f"[BLAST progress] '{self.file_path}' has {n_lines} lines so far...")
            else:
                logger.info("[BLAST progress] transcriptome.swissprot not created yet...")

    def stop(self):
        self._stop_event.set()


def run_gawn_with_monitor(gawn_command, file_path, use_conda, use_docker, output_dir):
    monitor = SwissprotMonitor(file_path=file_path, interval=60)
    monitor.start()
    try:
        run_pipeline_command(gawn_command, use_conda, use_docker, output_dir)
    finally:
        monitor.stop()
        monitor.join()

###############################################################################
# ENV FILES (UPDATED: includes portcullis, pasa, mikado; drops plotting libs/R)
###############################################################################
def environment_yml_content():
    return """\
name: annotate_env
channels:
  - bioconda
  - conda-forge
  - defaults
dependencies:
  # core
  - python=3.9
  - numpy
  - pandas
  - biopython
  - pybedtools
  - pysam
  - pyyaml
  - parallel
  - procps-ng
  - tqdm

  # IO / comparison
  - gffcompare
  - gffread
  - bedtools
  - samtools
  - seqkit

  # assembly/evidence
  - stringtie
  - blast
  - diamond
  - gmap

  # coding sequences
  - transdecoder

  # exon-correction toolchain
  - portcullis
  - mikado
  - pasa            # PASA pipeline with SQLite support
  - sqlite
"""

def dockerfile_content():
    return """\
FROM continuumio/miniconda3:4.8.2
WORKDIR /opt

# Create the environment
COPY environment.yml /tmp/environment.yml
RUN conda env create -f /tmp/environment.yml && conda clean -afy

# Make the env default
ENV PATH=/opt/conda/envs/annotate_env/bin:$PATH

# Sanity checks
RUN bash -lc "source activate annotate_env && which gffread && which portcullis && which mikado && which Launch_PASA_pipeline.pl"

WORKDIR /data
"""

def write_environment_yml():
    with open("environment.yml", "w") as f:
        f.write(environment_yml_content())
    log_green_info("environment.yml written.")

def write_dockerfile():
    with open("Dockerfile", "w") as f:
        f.write(dockerfile_content())
    log_green_info("Dockerfile written.")

###############################################################################
# Checking Tools
###############################################################################
def conda_run_which(tool):
    cmd = f"conda run -n annotate_env which {tool}"
    result = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return (result.returncode == 0)

def docker_run_which(tool):
    cmd = f"docker run --rm myorg/annotate_env:latest conda run -n annotate_env which {tool}"
    result = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return (result.returncode == 0)

def install_missing_tools_conda(missing):
    log_green_info("Installing missing tools in conda environment...")
    pkgs = " ".join(missing)
    cmd = f"conda run -n annotate_env conda install -c bioconda -c conda-forge -y {pkgs}"
    run_cmd(cmd)

def rebuild_docker_with_missing_tools(missing):
    log_green_info("Updating environment.yml to include missing tools for Docker...")
    with open("environment.yml","r") as f:
        lines=f.readlines()
    new_lines=[]
    insert_index=None
    for i,line in enumerate(lines):
        new_lines.append(line)
        if line.strip().startswith("dependencies:"):
            insert_index=i+1
    if insert_index is not None:
        for tool in missing:
            new_lines.insert(insert_index, f"  - {tool}\n")

    with open("environment.yml","w") as f:
        f.writelines(new_lines)
    log_green_info("Rebuilding Docker image with updated environment.yml...")
    write_dockerfile()
    run_cmd("docker build . -t myorg/annotate_env:latest")

def verify_tools_conda():
    log_green_info("Verifying tools in conda environment...")
    missing=[]
    for tool in REQUIRED_TOOLS:
        logger.info(f"Checking tool: {tool}")
        if not conda_run_which(tool):
            logger.warning(f"{tool} not found in conda environment.")
            missing.append(tool)
    if missing:
        install_missing_tools_conda(missing)
        for tool in missing:
            if not conda_run_which(tool):
                logger.error(f"{tool} still not found => exit.")
                sys.exit(1)
        log_green_info("All missing tools installed in conda.")
    else:
        log_green_info("All required tools present in conda environment.")

def verify_tools_docker():
    log_green_info("Verifying tools in docker image...")
    missing=[]
    for tool in REQUIRED_TOOLS:
        logger.info(f"Checking tool: {tool}")
        if not docker_run_which(tool):
            logger.warning(f"{tool} not found in Docker image.")
            missing.append(tool)
    if missing:
        log_green_info("Missing tools => rebuild docker environment.")
        rebuild_docker_with_missing_tools(missing)
        for tool in missing:
            if not docker_run_which(tool):
                logger.error(f"{tool} not found after rebuild => exit.")
                sys.exit(1)
        log_green_info("All missing tools present in Docker now.")
    else:
        log_green_info("All required tools present in Docker image.")

###############################################################################
# Env Installation
###############################################################################
def install_conda_env():
    log_green_info("Installing conda environment...")
    write_environment_yml()
    run_cmd("conda env create -f environment.yml")
    verify_tools_conda()
    log_green_info("::: Installation Complete. Exiting. :::")
    sys.exit(0)

def install_docker_image():
    log_green_info("Installing docker image...")
    write_environment_yml()
    write_dockerfile()
    run_cmd("docker build . -t myorg/annotate_env:latest")
    verify_tools_docker()
    log_green_info("::: Installation Complete. Exiting. :::")
    sys.exit(0)

def conda_env_exists():
    cmd="conda env list | grep annotate_env"
    r=subprocess.run(cmd, shell=True)
    return (r.returncode==0)

def docker_image_exists():
    cmd="docker images | grep myorg/annotate_env"
    r=subprocess.run(cmd, shell=True)
    return (r.returncode==0)

def purge_all_envs():
    logger.info("Purging conda env and docker image...")
    run_cmd("conda remove -n annotate_env --all -y || true")
    run_cmd("docker rmi myorg/annotate_env:latest -f || true")
    sys.exit(0)

###############################################################################
# Input handling
###############################################################################
def check_inputs(args):
    """
    Ensure output_dir is absolute, create if needed.
    Copy relevant input files to that folder so local references work:
      - egap_gff
      - a (GTF/GFF3)
      - g (FASTA)
      - bam
    Then reassign args.<input> to the local path.
    """
    if not args.output:
        logger.error("No --output directory provided; cannot proceed.")
        sys.exit(1)

    outdir = os.path.abspath(args.output)
    os.makedirs(outdir, exist_ok=True)
    args.output = outdir  # store updated

    if args.egap_gff and os.path.isfile(args.egap_gff):
        egap_abs = os.path.abspath(args.egap_gff)
        egap_local = os.path.join(outdir, os.path.basename(egap_abs))
        if not os.path.exists(egap_local):
            shutil.copy(egap_abs, egap_local)
        args.egap_gff = egap_local

    if args.a and os.path.isfile(args.a):
        a_abs = os.path.abspath(args.a)
        a_local = os.path.join(outdir, os.path.basename(a_abs))
        if not os.path.exists(a_local):
            shutil.copy(a_abs, a_local)
        args.a = a_local

    if args.g and os.path.isfile(args.g):
        g_abs = os.path.abspath(args.g)
        g_local = os.path.join(outdir, os.path.basename(g_abs))
        if not os.path.exists(g_local):
            shutil.copy(g_abs, g_local)
        args.g = g_local

    if args.bam and os.path.isfile(args.bam):
        bam_abs = os.path.abspath(args.bam)
        bam_local = os.path.join(outdir, os.path.basename(bam_abs))
        if not os.path.exists(bam_local):
            shutil.copy(bam_abs, bam_local)
        args.bam = bam_local

###############################################################################
# PREPROCESSING (kept): GTF normalization + optional name merge with EGAP GFF
###############################################################################
def run_preprocessing(args):
    """
    Minimal preprocessing (kept):
      * If --egap_gff not provided:
          - run unique_gene_id.py on -a → preprocessed_best.gtf
      * If --egap_gff provided:
          - run unique_gene_id.py on -a
          - merge names with merge_stringtie_names.py using --egap_gff
          - output => preprocessed_best.gtf
    No StringTie parameter tuning is performed anymore.
    """
    outdir = args.output

    if not args.a or not os.path.isfile(args.a):
        logger.error("[PREPROCESSING] requires -a <some.gtf/gff3>")
        sys.exit(1)

    # Fetch helpers
    cmd_wget_uniq = (
        "wget -O unique_gene_id.py "
        "https://raw.githubusercontent.com/cfarkas/amyg/refs/heads/main/third_parties/unique_gene_id.py"
    )
    run_pipeline_command(cmd_wget_uniq, args.use_conda, args.use_docker, outdir)
    run_pipeline_command("chmod 755 unique_gene_id.py", args.use_conda, args.use_docker, outdir)

    a_basename = os.path.basename(args.a)
    a_base_noext = os.path.splitext(a_basename)[0]
    uniq_out = a_base_noext + ".unique_gene_id.gtf"

    # Convert to GTF if needed, then assign unique gene ids
    tmp_gtf = a_basename
    if a_basename.lower().endswith(".gff3") or a_basename.lower().endswith(".gff"):
        # normalize to GTF for downstream helpers
        tmp_gtf = a_base_noext + ".as_gtf.gtf"
        run_pipeline_command(f"gffread -T -o {tmp_gtf} {a_basename}", args.use_conda, args.use_docker, outdir)

    run_pipeline_command(f"python unique_gene_id.py {tmp_gtf}", args.use_conda, args.use_docker, outdir)
    if not os.path.isfile(os.path.join(outdir, uniq_out)):
        # try fallback glob
        alt = glob.glob(os.path.join(outdir, a_base_noext + ".unique_gene_id.gtf"))
        if not alt:
            logger.error("[PREPROCESSING] unique_gene_id failed.")
            sys.exit(1)

    final_pre = os.path.join(outdir, "preprocessed_best.gtf")

    if args.egap_gff:
        logger.info("=== Preprocessing: Downloading merge_stringtie_names.py ===")
        run_pipeline_command(
            "wget https://raw.githubusercontent.com/cfarkas/amyg/refs/heads/main/scripts/merge_stringtie_names.py",
            args.use_conda,
            args.use_docker,
            outdir
        )
        run_pipeline_command("chmod 755 merge_stringtie_names.py", args.use_conda, args.use_docker, outdir)
        out_named = "transcripts_named.gtf"
        cmd_merge = (
            f"python merge_stringtie_names.py "
            f"--stringtie_gtf {uniq_out} "
            f"--egap_gff {os.path.basename(args.egap_gff)} "
            f"--output_gtf {out_named}"
        )
        run_pipeline_command(cmd_merge, args.use_conda, args.use_docker, outdir)
        shutil.copy(os.path.join(outdir, out_named), final_pre)
    else:
        shutil.copy(os.path.join(outdir, uniq_out), final_pre)

    logger.info(f"[PREPROCESSING] => Overriding -a => {final_pre}")
    args.a = final_pre

###############################################################################
# EXON CORRECTION (NEW): Portcullis + (PASA | Mikado)
###############################################################################
def find_portcullis_pass_bed(pc_dir):
    """Return path to filtered/pass junctions bed from Portcullis."""
    candidates = [
        os.path.join(pc_dir, "filtered_junctions.bed"),
        os.path.join(pc_dir, "pass.junctions.bed"),
        os.path.join(pc_dir, "junc/filtered_junctions.bed"),
    ]
    for p in candidates:
        if os.path.isfile(p):
            return p
    return None

def write_pasa_align_config(outdir, use_docker):
    """
    Minimal PASA SQLite config.
    Creates: pasa.alignAssembly.config in outdir, using local path 'pasa.sqlite'
    (absolute path inside container is /data/pasa.sqlite).
    """
    conf_path = os.path.join(outdir, "pasa.alignAssembly.config")
    lines = [
        "# PASA SQLite config (auto-generated by amyg.py)",
        "DATABASE=pasa.sqlite",
        "",
        "# Validation thresholds:",
        "validate_alignments_in_db.dbi:--MIN_PERCENT_ALIGNED=80",
        "validate_alignments_in_db.dbi:--MIN_AVG_PER_ID=80",
    ]
    with open(conf_path, "w") as f:
        f.write("\n".join(lines) + "\n")
    return conf_path

def run_portcullis(args):
    outdir = args.output
    g = os.path.basename(args.g)
    b = os.path.basename(args.bam)
    pc_dir = os.path.join(outdir, "portcullis_out")
    os.makedirs(pc_dir, exist_ok=True)
    cmd = f"portcullis full -o {pc_dir} -t {args.threads} {g} {b}"
    log_green_info(f"[EXON_CORR] Portcullis: {cmd}")
    run_pipeline_command(cmd, args.use_conda, args.use_docker, outdir)
    bed = find_portcullis_pass_bed(pc_dir)
    if not bed:
        logger.error("[EXON_CORR] Portcullis pass junctions BED not found.")
        sys.exit(1)
    return bed, pc_dir

def run_exon_correction_pasa(args, portcullis_bed):
    """
    PASA alignment assembly to refine exon boundaries based on transcript alignments.
    Input transcripts are generated from current args.a using gffread.
    """
    outdir = args.output
    g = os.path.basename(args.g)
    a = os.path.basename(args.a)

    # Generate transcripts.fa from current annotation
    run_pipeline_command(
        f"gffread -w transcripts_for_pasa.fa -g {g} {a}",
        args.use_conda, args.use_docker, outdir
    )

    conf = write_pasa_align_config(outdir, args.use_docker)
    pasa_log = os.path.join(outdir, "pasa.align.log")
    # Run PASA (GMAP-based)
    cmd = (
        f"Launch_PASA_pipeline.pl "
        f"-c {os.path.basename(conf)} -C -R "
        f"-g {g} -t transcripts_for_pasa.fa "
        f"--ALIGNERS gmap "
        f"--CPU {args.threads} "
        f"> {os.path.basename(pasa_log)} 2>&1"
    )
    log_green_info("[EXON_CORR] Running PASA alignment assembly...")
    run_pipeline_command(cmd, args.use_conda, args.use_docker, outdir)

    # Detect PASA GFF3 output
    candidates = sorted(
        glob.glob(os.path.join(outdir, "*.pasa_assemblies.gff3"))
        + glob.glob(os.path.join(outdir, "*.gene_structures_post_PASA_updates*.gff3"))
    )
    if not candidates:
        logger.error("[EXON_CORR] PASA GFF3 not found. Inspect pasa.align.log.")
        sys.exit(1)

    pasa_gff3 = os.path.basename(candidates[-1])
    corrected_gtf = "pasa_corrected.gtf"
    run_pipeline_command(
        f"gffread -T -o {corrected_gtf} {pasa_gff3}",
        args.use_conda, args.use_docker, outdir
    )
    logger.info(f"[EXON_CORR] PASA produced {pasa_gff3} -> {corrected_gtf}")
    args.a = os.path.join(outdir, corrected_gtf)

def run_exon_correction_mikado(args, portcullis_bed):
    """
    Mikado: configure → prepare → (optional) serialise → pick
    Uses current args.a (GTF/GFF3) and genome, integrates Portcullis junctions.
    """
    outdir = args.output
    g = os.path.basename(args.g)
    a = os.path.basename(args.a)

    mk_dir = os.path.join(outdir, "mikado")
    os.makedirs(mk_dir, exist_ok=True)

    config_yaml = os.path.join("mikado", "mikado_config.yaml")
    # configure
    cmd_cfg = (
        f"mikado configure "
        f"--gff {a} "
        f"--reference {g} "
        f"--junctions {os.path.relpath(portcullis_bed, outdir)} "
        f"-y -od mikado mikado_config.yaml"
    )
    log_green_info("[EXON_CORR] Mikado configure")
    run_pipeline_command(cmd_cfg, args.use_conda, args.use_docker, outdir)

    # prepare
    cmd_prep = (
        f"mikado prepare "
        f"--json-conf {config_yaml} "
        f"--reference {g} "
        f"-od mikado "
        f"-o mikado_prepared.gtf -of mikado_prepared.fasta "
        f"{a}"
    )
    run_pipeline_command(cmd_prep, args.use_conda, args.use_docker, outdir)

    # optional: ORFs (can improve pick, but not strictly required)
    if os.path.isfile(os.path.join(outdir, "mikado", "mikado_prepared.fasta")):
        run_pipeline_command(
            "TransDecoder.LongOrfs -t mikado/mikado_prepared.fasta",
            args.use_conda, args.use_docker, outdir
        )
        run_pipeline_command(
            "TransDecoder.Predict -t mikado/mikado_prepared.fasta",
            args.use_conda, args.use_docker, outdir
        )

    # serialise (will load junctions from config; ORF/homology if available)
    cmd_ser = f"mikado serialise --json-conf {config_yaml}"
    run_pipeline_command(cmd_ser, args.use_conda, args.use_docker, outdir)

    # pick
    cmd_pick = f"mikado pick --json-conf {config_yaml} -o mikado/mikado.loci.gff3"
    run_pipeline_command(cmd_pick, args.use_conda, args.use_docker, outdir)

    loci_gff3 = os.path.join(outdir, "mikado", "mikado.loci.gff3")
    if not os.path.isfile(loci_gff3):
        logger.error("[EXON_CORR] Mikado did not produce mikado.loci.gff3")
        sys.exit(1)

    corrected_gtf = os.path.join(outdir, "mikado_corrected.gtf")
    run_pipeline_command(
        f"gffread -T -o {os.path.basename(corrected_gtf)} mikado/mikado.loci.gff3",
        args.use_conda, args.use_docker, outdir
    )
    logger.info(f"[EXON_CORR] Mikado produced mikado.loci.gff3 -> {corrected_gtf}")
    args.a = corrected_gtf

def run_exon_correction(args):
    """
    Orchestrate exon-correction if requested.
    Requires: --bam and -g when --exon_correction != none
    """
    if args.exon_correction == "none":
        logger.info("[EXON_CORR] Skipping exon correction (none).")
        return

    if not args.bam or not os.path.isfile(args.bam):
        logger.error("[EXON_CORR] --bam is required for Portcullis.")
        sys.exit(1)
    if not args.g or not os.path.isfile(args.g):
        logger.error("[EXON_CORR] -g genome FASTA is required.")
        sys.exit(1)
    if not args.a or not os.path.isfile(args.a):
        logger.error("[EXON_CORR] -a annotation (GTF/GFF3) is required.")
        sys.exit(1)

    # Always run Portcullis first
    junc_bed, pc_dir = run_portcullis(args)

    # Then run the chosen corrector
    if args.exon_correction == "pasa":
        run_exon_correction_pasa(args, junc_bed)
    elif args.exon_correction == "mikado":
        run_exon_correction_mikado(args, junc_bed)
    else:
        logger.error(f"[EXON_CORR] Unknown option: {args.exon_correction}")
        sys.exit(1)

###############################################################################
# SINGLE-CELL LOGIC (unchanged)
###############################################################################
def run_single_cell_mode(args):
    """
    Single-cell pipeline that:
      1) Copies .bam + ref GTF/FA => outdir.
      2) For each .bam => run StringTie => produce GTF.
      3) Merge all => all_merged.gtf (stringtie --merge).
      4) unique_gene_id.py + merge_stringtie_names.py => final named GTF.
      5) Filter novel => transcripts_truly_novel.gtf.
      6) Update args.a => that novel GTF, args.g => local reference FA, args.output => single-cell folder.
      7) Return so the normal pipeline can continue inline with valid -a, -g, -o.
    """
    logger.info("[SINGLE_CELL] Starting single-cell pipeline...")

    date_str = datetime.datetime.now().strftime("%Y%m%d")
    outdir = os.path.join(args.output, f"amyg_singlecell_{date_str}")
    os.makedirs(outdir, exist_ok=True)

    # Basic checks
    if not args.input_dir or not args.ref_gtf or not args.ref_fa:
        logger.error("[SINGLE_CELL] => requires --input_dir, --ref_gtf, --ref_fa")
        sys.exit(1)

    # 1) Copy reference GTF/FA => outdir
    ref_gtf_abs = os.path.abspath(args.ref_gtf)
    ref_gtf_basename = os.path.basename(ref_gtf_abs)
    local_ref_gtf = os.path.join(outdir, ref_gtf_basename)
    if not os.path.exists(local_ref_gtf):
        shutil.copy(ref_gtf_abs, local_ref_gtf)

    ref_fa_abs = os.path.abspath(args.ref_fa)
    ref_fa_basename = os.path.basename(ref_fa_abs)
    local_ref_fa = os.path.join(outdir, ref_fa_basename)
    if not os.path.exists(local_ref_fa):
        shutil.copy(ref_fa_abs, local_ref_fa)

    # 2) Copy .bam => outdir => run StringTie => single GTF
    bam_list = glob.glob(os.path.join(os.path.abspath(args.input_dir), "*.bam"))
    if not bam_list:
        logger.error("[SINGLE_CELL] => no .bam found => abort")
        sys.exit(1)

    local_bams = []
    for bam_file in bam_list:
        base = os.path.basename(bam_file)
        local_bam = os.path.join(outdir, base)
        if not os.path.exists(local_bam):
            shutil.copy(bam_file, local_bam)
        local_bams.append(local_bam)

    single_gtfs = []
    for local_bam in local_bams:
        bam_name = os.path.basename(local_bam)
        out_gtf = os.path.join(outdir, bam_name + ".stringtie.gtf")
        cmd_st = [
            "stringtie",
            "-p", str(args.threads),
            "-G", ref_gtf_basename,  # local ref GTF
            "-c","2","-s","8","-f","0.05","-j","3","-m","300","-v",
            "-o", os.path.basename(out_gtf),
            bam_name
        ]
        run_pipeline_command(" ".join(cmd_st), args.use_conda, args.use_docker, outdir)
        if os.path.isfile(out_gtf):
            single_gtfs.append(out_gtf)
        else:
            logger.warning(f"[SINGLE_CELL] Missing {out_gtf}")

    if not single_gtfs:
        logger.warning("[SINGLE_CELL] => no GTF => skip pipeline.")
        sys.exit(0)

    # 3) stringtie --merge => all_merged.gtf
    all_merged = os.path.join(outdir, "all_merged.gtf")
    list_file = os.path.join(outdir, "gtf_list.txt")
    with open(list_file, "w") as lf:
        for gtf_path in single_gtfs:
            lf.write(os.path.basename(gtf_path) + "\n")

    merge_cmd = [
        "stringtie", "--merge",
        "-l", "STRG",
        "-G", ref_gtf_basename,
        "-o", "all_merged.gtf",
        os.path.basename(list_file)
    ]
    run_pipeline_command(" ".join(merge_cmd), args.use_conda, args.use_docker, outdir)

    if not os.path.isfile(all_merged):
        logger.error("[SINGLE_CELL] => missing all_merged.gtf => merge failed.")
        sys.exit(1)

    # 4) unique_gene_id.py + merge_stringtie_names.py => final named GTF
    unique_py = os.path.join(outdir, "unique_gene_id.py")
    if not os.path.isfile(unique_py):
        cmd_wget_uniq = (
            "wget -O unique_gene_id.py "
            "https://raw.githubusercontent.com/cfarkas/amyg/refs/heads/main/third_parties/unique_gene_id.py"
        )
        run_pipeline_command(cmd_wget_uniq, args.use_conda, args.use_docker, outdir)
        run_pipeline_command("chmod 755 unique_gene_id.py", args.use_conda, args.use_docker, outdir)

    merge_py = os.path.join(outdir, "merge_stringtie_names.py")
    if not os.path.isfile(merge_py):
        cmd_wget_merge = (
            "wget https://raw.githubusercontent.com/cfarkas/amyg/refs/heads/main/scripts/merge_stringtie_names.py"
        )
        run_pipeline_command(cmd_wget_merge, args.use_conda, args.use_docker, outdir)
        run_pipeline_command("chmod 755 merge_stringtie_names.py", args.use_conda, args.use_docker, outdir)

    base_noext = os.path.splitext(all_merged)[0]
    uniq_gtf   = base_noext + ".unique_gene_id.gtf"
    cmd_uniq   = f"python {os.path.basename(unique_py)} all_merged.gtf"
    run_pipeline_command(cmd_uniq, args.use_conda, args.use_docker, outdir)

    if not os.path.isfile(uniq_gtf):
        logger.warning("[SINGLE_CELL] => unique_gene_id => missing => skip naming.")
        final_named = all_merged
    else:
        named_gtf = base_noext + "_named.gtf"
        cmd_merge_names = (
            f"python {os.path.basename(merge_py)} "
            f"--stringtie_gtf {os.path.basename(uniq_gtf)} "
            f"--egap_gff {ref_gtf_basename} "
            f"--output_gtf {os.path.basename(named_gtf)}"
        )
        run_pipeline_command(cmd_merge_names, args.use_conda, args.use_docker, outdir)
        final_named = named_gtf if os.path.isfile(named_gtf) else uniq_gtf

    # 5) Filter novel => transcripts_truly_novel.gtf
    def parse_attrs_local(a_str):
        dd = {}
        for chunk in a_str.split(';'):
            chunk = chunk.strip()
            if not chunk:
                continue
            parts = chunk.split(' ', 1)
            if len(parts) < 2:
                continue
            k, v = parts
            v = v.strip().strip('"')
            dd[k] = v
        return dd

    truly_novel = os.path.join(outdir, "transcripts_truly_novel.gtf")

    def produce_truly_novel(gtf_in, novel_out):
        if not os.path.isfile(gtf_in):
            logger.warning(f"[SINGLE_CELL] produce_truly_novel => missing {gtf_in}")
            return
        lines_in, lines_novel = 0, 0
        with open(gtf_in, "r") as fin, open(novel_out, "w") as fout:
            for line in fin:
                if line.startswith("#") or not line.strip():
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 9:
                    continue
                attrs = parse_attrs_local(fields[8])
                gene_id = attrs.get("gene_id", "")
                lines_in += 1
                if gene_id.startswith("STRG."):
                    lines_novel += 1
                    fout.write(line)
        logger.info(
            f"[SINGLE_CELL] produce_truly_novel => scanned {lines_in}, "
            f"found {lines_novel} novel => {novel_out}"
        )

    produce_truly_novel(final_named, truly_novel)

    # 6) Reassign args so the main pipeline sees -a, -g, and -o
    logger.info("[SINGLE_CELL] => Done producing transcripts_truly_novel.gtf.")
    args.a = truly_novel
    args.g = local_ref_fa
    args.output = outdir

    logger.info("[SINGLE_CELL] => returning to main => normal pipeline will continue with -a, -g, -o set properly.")
    return

###############################################################################
# MAIN
###############################################################################
def main():
    parser = argparse.ArgumentParser(
        description="annotation pipeline with optional single_cell and exon-correction (Portcullis + PASA/Mikado)."
    )
    parser.add_argument("--install", choices=["conda", "docker"], help="Install environment and exit.")
    parser.add_argument("--use_conda", action="store_true", help="Use conda env.")
    parser.add_argument("--use_docker", action="store_true", help="Use docker image.")
    parser.add_argument("--threads", type=int, default=10)
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--purge_all_envs", action="store_true")

    parser.add_argument("-o", "--output", help="Output directory")
    parser.add_argument("-a", help="Input annotation (GTF/GFF3)")
    parser.add_argument("-g", help="Reference genome (FASTA)")
    parser.add_argument("--egap_gff", help="EGAP GFF for merging or reference checks")

    # single-cell mode (kept)
    parser.add_argument("--single_cell", action="store_true",
                        help="If set => multi-bam single-cell logic => exit after done.")
    parser.add_argument("--input_dir", help="Directory with .bam for single_cell")
    parser.add_argument("--ref_gtf", help="Known ref GTF for single_cell mode")
    parser.add_argument("--ref_fa", help="Reference genome for single_cell final pass")
    parser.add_argument("--overlap_frac", type=float, default=0.05,
                        help="Overlap fraction for single_cell mode.")

    # simple preprocessing (kept; no tuning)
    parser.add_argument("--preprocessing", action="store_true",
                        help="Normalize IDs and optionally merge names with --egap_gff (no StringTie tuning).")

    # exon-correction (NEW)
    parser.add_argument("--exon_correction", choices=["none", "pasa", "mikado"], default="none",
                        help="Run Portcullis plus PASA or Mikado to refine exon boundaries (requires --bam and -g).")
    parser.add_argument("--bam", help="RNA-seq BAM for Portcullis (required when --exon_correction != none)")

    args = parser.parse_args()

    if args.purge_all_envs:
        purge_all_envs()
    if args.install == "conda":
        install_conda_env()
    elif args.install == "docker":
        install_docker_image()

    # unify absolute paths and local copies
    check_inputs(args)

    # single-cell branch (produces args.a/args.g and returns)
    if args.single_cell:
        run_single_cell_mode(args)

    # light preprocessing (kept)
    if args.preprocessing:
        run_preprocessing(args)

    # optional exon-correction (Portcullis + PASA/Mikado)
    if args.exon_correction != "none":
        run_exon_correction(args)

    # === Normal pipeline continues ===
    log_green_info("Starting normal pipeline script...")

    if not args.a or not args.g or not args.output:
        logger.error("For normal pipeline => must provide -a, -g, -o.")
        sys.exit(1)

    use_conda = args.use_conda
    use_docker = args.use_docker
    threads = args.threads
    force = args.force
    a = args.a
    g = args.g
    output_dir = args.output

    # Check environment usage
    if use_conda and not conda_env_exists():
        logger.error("Conda env 'annotate_env' not found => run --install conda first.")
        sys.exit(1)
    if use_docker and not docker_image_exists():
        logger.error("Docker image 'myorg/annotate_env:latest' not found => run --install docker first.")
        sys.exit(1)

    if not os.path.isdir(output_dir):
        logger.error(f"Output dir => {output_dir} not found => create it first.")
        sys.exit(1)

    db_dir = os.path.join(output_dir, "database")
    gawn_config_path = os.path.join(output_dir, "gawn_config.sh")
    if (os.path.exists(db_dir) or os.path.exists(gawn_config_path)) and not force:
        logger.error("db or gawn_config.sh exist => use --force to overwrite.")
        sys.exit(1)
    if force:
        if os.path.exists(db_dir):
            shutil.rmtree(db_dir)
        if os.path.exists(gawn_config_path):
            os.remove(gawn_config_path)

    log_green_info("Downloading & preparing SwissProt inside output_dir/database...")
    os.makedirs(db_dir, exist_ok=True)
    run_cmd(f"wget -P {db_dir} ftp://ftp.ncbi.nlm.nih.gov/blast/db/swissprot.tar.gz")
    run_cmd(f"gunzip {os.path.join(db_dir,'swissprot.tar.gz')}")
    run_cmd(f"tar -xvf {os.path.join(db_dir,'swissprot.tar')} -C {db_dir}")

    with open(gawn_config_path, "w") as gc:
        gc.write("#!/bin/bash\n")
        gc.write(f"NCPUS={threads}\n")
        gc.write("SKIP_GENOME_INDEXING=1\n")
        gc.write("GENOME_NAME=\"genome.fasta\"\n")
        gc.write("TRANSCRIPTOME_NAME=\"transcriptome.fasta\"\n")
        if use_docker:
            gc.write('SWISSPROT_DB="/data/database/swissprot"\n')
        else:
            gc.write('SWISSPROT_DB="03_data/swissprot"\n')
        gc.write("#\n")

    a_abs = os.path.abspath(a)
    g_abs = os.path.abspath(g)
    a_filename = os.path.basename(a_abs)
    g_filename = os.path.basename(g_abs)

    logger.info("::: Step 1 => gffread => transcripts.fa :::")
    run_pipeline_command(
        f"gffread -w transcripts.fa -g {g_filename} {a_filename}",
        use_conda, use_docker, output_dir
    )

    logger.info("::: Step 2 => GAWN for gene annotation :::")
    GAWN_DIR = os.path.join(output_dir, "gawn")
    if os.path.isdir(GAWN_DIR):
        shutil.rmtree(GAWN_DIR)
    os.makedirs(GAWN_DIR, exist_ok=True)
    run_pipeline_command("git clone https://github.com/enormandeau/gawn.git gawn",
                         use_conda, use_docker, output_dir)

    shutil.copy(os.path.join(output_dir, "transcripts.fa"),
                os.path.join(GAWN_DIR, "03_data", "transcriptome.fasta"))
    shutil.copy(os.path.join(output_dir, g_filename),
                os.path.join(GAWN_DIR, "03_data", "genome.fasta"))

    if use_conda:
        def copy_swissprot_conda(db_dir_, gawn_dir_):
            logger.info("Copying SwissProt => gawn/03_data => conda mode")
            three_data = os.path.join(gawn_dir_, "03_data")
            os.makedirs(three_data, exist_ok=True)
            run_cmd(f"cp -v {db_dir_}/swissprot.* {three_data}/")

        copy_swissprot_conda(db_dir, GAWN_DIR)
        with open(gawn_config_path, "r") as f:
            lines = f.readlines()
        new_lines = []
        for line in lines:
            if "SWISSPROT_DB" in line:
                new_lines.append('SWISSPROT_DB="03_data/swissprot"\n')
            else:
                new_lines.append(line)
        with open(gawn_config_path, "w") as f:
            f.writelines(new_lines)

    shutil.copy(gawn_config_path, os.path.join(GAWN_DIR, "02_infos", "gawn_config.sh"))
    swissprot_file = os.path.join(GAWN_DIR, "04_annotation", "transcriptome.swissprot")
    logger.info("::: Running GAWN pipeline + progress monitor :::")
    run_gawn_with_monitor("cd gawn && ./gawn 02_infos/gawn_config.sh",
                          swissprot_file,
                          use_conda, use_docker,
                          output_dir)

    hits_path = os.path.join(GAWN_DIR, "04_annotation", "transcriptome.hits")
    swissprot_path = os.path.join(GAWN_DIR, "04_annotation", "transcriptome.swissprot")
    if not os.path.isfile(hits_path) or not os.path.isfile(swissprot_path):
        logger.error("transcriptome.hits or transcriptome.swissprot not found => GAWN error.")
        sys.exit(9999)
    annotation_table_path = os.path.join(GAWN_DIR, "05_results", "transcriptome_annotation_table.tsv")
    if not os.path.isfile(annotation_table_path):
        logger.error("transcriptome_annotation_table.tsv missing => GAWN error.")
        sys.exit(9999)
    shutil.copy(swissprot_path, output_dir)
    shutil.copy(hits_path, output_dir)
    shutil.copy(annotation_table_path, output_dir)

    logger.info("::: Step 4 => TransDecoder :::")
    td_dir = os.path.join(output_dir, "transcripts.fa.transdecoder_dir")
    if os.path.exists(td_dir):
        shutil.rmtree(td_dir)

    run_pipeline_command("TransDecoder.LongOrfs -t transcripts.fa",
                         use_conda, use_docker, output_dir)
    run_pipeline_command("TransDecoder.Predict -t transcripts.fa",
                         use_conda, use_docker, output_dir)

    logger.info("::: Step 5 => final_results :::")
    FINAL_RESULTS_DIR = os.path.join(output_dir, "final_results")
    os.makedirs(FINAL_RESULTS_DIR, exist_ok=True)
    to_final = [
        "transcripts.fa.transdecoder_dir/longest_orfs.gff3",
        "transcripts.fa.transdecoder_dir/longest_orfs.cds",
        "transcripts.fa.transdecoder_dir/longest_orfs.pep",
        "transcriptome.hits",
        "transcriptome.swissprot",
        "transcripts.fa",
        "transcriptome_annotation_table.tsv"
    ]
    for f in to_final:
        src = os.path.join(output_dir, f)
        if os.path.isfile(src):
            shutil.move(src, FINAL_RESULTS_DIR)

    logger.info("::: Step 6 => transdecoder_results :::")
    TRANSDECODER_RESULTS_DIR = os.path.join(output_dir, "transdecoder_results")
    os.makedirs(TRANSDECODER_RESULTS_DIR, exist_ok=True)

    for f in os.listdir(output_dir):
        if f.startswith("transcripts.fa.transdecoder."):
            shutil.move(os.path.join(output_dir, f), TRANSDECODER_RESULTS_DIR)

    for f in os.listdir(os.getcwd()):
        if f.startswith("pipeliner.") and f.endswith(".cmds"):
            shutil.move(os.path.join(os.getcwd(), f), TRANSDECODER_RESULTS_DIR)

    leftover_dirs = [
        "transcripts.fa.transdecoder_dir",
        "transcripts.fa.transdecoder_dir.__checkpoints",
        "transcripts.fa.transdecoder_dir.__checkpoints_longorfs"
    ]
    for d in leftover_dirs:
        dpath = os.path.join(os.getcwd(), d)
        if os.path.isdir(dpath):
            shutil.move(dpath, TRANSDECODER_RESULTS_DIR)
    logger.info("::: TransDecoder => transdecoder_results dir :::")

    logger.info("::: Step 7 => annotate GTF :::")
    logger.info("::: Downloading annotate_gtf.py :::")
    run_pipeline_command(
        "curl -O https://raw.githubusercontent.com/cfarkas/amyg/refs/heads/main/scripts/annotate_gtf.py",
        use_conda, use_docker, output_dir
    )
    run_pipeline_command("chmod +x annotate_gtf.py", use_conda, use_docker, output_dir)

    gtf_basename = os.path.basename(a_abs)
    INPUT_GTF = gtf_basename
    HITS_FILE = "gawn/04_annotation/transcriptome.hits"
    ANNOTATION_TABLE = "gawn/05_results/transcriptome_annotation_table.tsv"
    TEMP_ANNOTATED_GTF = "final_annotated.gtf"

    host_gtf_path = os.path.join(output_dir, gtf_basename)
    if not os.path.isfile(host_gtf_path):
        logger.error(f"Input GTF not found => {host_gtf_path}")
        sys.exit(9999)
    host_hits_path = os.path.join(output_dir, HITS_FILE)
    if not os.path.isfile(host_hits_path):
        logger.error(f"Hits file not found => {host_hits_path}")
        sys.exit(9999)
    host_table_path = os.path.join(output_dir, ANNOTATION_TABLE)
    if not os.path.isfile(host_table_path):
        logger.error(f"annotation_table not found => {host_table_path}")
        sys.exit(9999)

    run_pipeline_command(
        f"python annotate_gtf.py {INPUT_GTF} {HITS_FILE} {ANNOTATION_TABLE} {TEMP_ANNOTATED_GTF}",
        use_conda, use_docker, output_dir
    )
    logger.info("::: GTF Annotation Completed :::")

    local_annotated_gtf = os.path.join(output_dir, "final_annotated.gtf")
    if os.path.isfile(local_annotated_gtf):
        shutil.move(local_annotated_gtf, FINAL_RESULTS_DIR)
        logger.info("::: final_annotated.gtf => final_results :::")
    else:
        logger.warning("No final_annotated.gtf => check annotate_gtf.py logs")

    # wrap up directory layout
    contents = os.listdir(output_dir)
    exclude = {
        "final_results",
        "transdecoder_results",
        "database",
        "gawn_config.sh",
        "transcripts.fa",
        "gawn",
        "mikado",
        "portcullis_out",
        "pasa.alignAssembly.config",
        "pasa.sqlite",
        "pasa.align.log",
    }
    leftover = [c for c in contents if c not in exclude]
    if leftover:
        TIMESTAMP = time.strftime("%Y%m%d_%H%M%S")
        NEW_DIR = os.path.join(output_dir, f"amyg_{TIMESTAMP}")
        os.makedirs(NEW_DIR, exist_ok=True)
        logger.info(f"moving leftover => {NEW_DIR}")
        if os.path.exists(os.path.join(output_dir, "final_results")):
            shutil.move(os.path.join(output_dir, "final_results"), os.path.join(NEW_DIR, "final_results"))
        if os.path.exists(os.path.join(output_dir, "transdecoder_results")):
            shutil.move(os.path.join(output_dir, "transdecoder_results"), os.path.join(NEW_DIR, "transdecoder_results"))
        FINAL_DIR = os.path.join(NEW_DIR, "final_results")
    else:
        FINAL_DIR = os.path.join(output_dir, "final_results")

    logger.info("::: Pipeline completed :::")


if __name__ == "__main__":
    main()
