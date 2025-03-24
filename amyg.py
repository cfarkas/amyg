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
import glob
import statistics
import datetime

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
            f"bash -c \"{cmd}\""
        )
    elif use_conda:
        full_cmd = f"conda run -n annotate_env bash -c \"cd {output_dir} && {cmd}\""
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
# ENV FILES
###############################################################################
def environment_yml_content(use_docker=False):
    if not use_docker:
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
  - pybedtools
"""
    else:
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
  - pybedtools
"""

def dockerfile_content(use_docker=False):
    if not use_docker:
        return """\
FROM continuumio/miniconda3:4.8.2
WORKDIR /opt

COPY environment.yml /tmp/environment.yml
RUN conda env create -f /tmp/environment.yml && conda clean -afy

ENV PATH /opt/conda/envs/annotate_env/bin:$PATH
WORKDIR /data
"""
    else:
        return """\
FROM continuumio/miniconda3:4.8.2
WORKDIR /opt

COPY environment.yml /tmp/environment.yml

# 1) Create the conda environment
RUN conda env create -f /tmp/environment.yml && conda clean -afy

# 2) Install compilers in the environment so pip can build C/C++ packages
RUN conda run -n annotate_env conda install -c conda-forge -c bioconda -y compilers

# 3) Pip install Python data libs
RUN conda run -n annotate_env pip install numpy pandas matplotlib seaborn biopython intervaltree

# 4) Basic check
RUN conda run -n annotate_env python -c "import numpy, pandas, matplotlib, seaborn, Bio, intervaltree; print('Python packages installed correctly.')"

ENV PATH /opt/conda/envs/annotate_env/bin:$PATH
WORKDIR /data
"""

def write_environment_yml(use_docker):
    env_content = environment_yml_content(use_docker=use_docker)
    with open("environment.yml", "w") as f:
        f.write(env_content)
    log_green_info(f"environment.yml written (use_docker={use_docker}).")

def write_dockerfile(use_docker, extra_files=None):
    df_content = dockerfile_content(use_docker=use_docker).splitlines()
    with open("Dockerfile", "w") as f:
        f.write("\n".join(df_content) + "\n")
    log_green_info(f"Dockerfile written (use_docker={use_docker}).")

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
    write_dockerfile(use_docker=True)
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
    write_environment_yml(use_docker=False)
    run_cmd("conda env create -f environment.yml")
    verify_tools_conda()
    log_green_info("::: Installation Complete. Exiting. :::")
    sys.exit(0)

def install_docker_image():
    log_green_info("Installing docker image...")
    write_environment_yml(use_docker=True)
    write_dockerfile(use_docker=True)
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
# check_inputs
###############################################################################
def check_inputs(args):
    """
    Ensure output_dir is absolute, create if needed.
    Copy relevant input files to that folder so local references work:
      - egap_gff
      - a (GTF)
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

    # If we have a known file, copy it into outdir if not existing
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
# F1 from GFFcompare
###############################################################################
def compute_f1_from_stats(stats_file):
    if not os.path.isfile(stats_file):
        logger.warning(f"compute_f1_from_stats => missing {stats_file}")
        return 0.0

    sensitivity = 0.0
    precision   = 0.0
    found = False
    with open(stats_file) as sf:
        for line in sf:
            line=line.strip()
            if line.startswith("Transcript level:"):
                parts=line.split()
                if len(parts) >= 5:
                    try:
                        sens = float(parts[2])
                        prec = float(parts[4])
                        sensitivity = sens
                        precision   = prec
                        found = True
                        break
                    except:
                        pass

    if not found:
        return 0.0
    if (sensitivity+precision)==0:
        return 0.0
    f1 = 2.0*(sensitivity*precision)/(sensitivity+precision)
    logger.info(f"[F1] transcript-level => S={sensitivity:.2f}, P={precision:.2f}, F1={f1:.2f}")
    return f1

###############################################################################
# detect_coverage for bam => average coverage
###############################################################################
def detect_coverage(bam, outdir, use_conda, use_docker):
    logger.info("[PREPROCESSING] Detect coverage with samtools depth -a -d 0 ...")
    coverage_file = "temp_coverage.txt"
    bam_abs = os.path.abspath(bam)
    cmd_depth = f"samtools depth -a -d 0 {bam_abs} > {coverage_file}"
    run_pipeline_command(cmd_depth, use_conda, use_docker, outdir)

    full_cov = os.path.join(outdir, coverage_file)
    if not os.path.isfile(full_cov):
        logger.warning("[COVERAGE] coverage file not found => return None")
        return None

    covs=[]
    with open(full_cov) as f:
        for line in f:
            parts=line.strip().split()
            if len(parts)==3:
                try:
                    c=float(parts[2])
                    covs.append(c)
                except:
                    pass

    if not covs:
        logger.warning("[COVERAGE] no coverage lines => None")
        return None

    avg=statistics.mean(covs)
    logger.info(f"[COVERAGE] average coverage => {avg:.2f}")
    return avg


###############################################################################
# Preprocessing
###############################################################################

def run_preprocessing_tuning(args):
    import glob
    import os
    import sys
    import shutil
    import logging

    outdir = args.output  # assume already an absolute path

    # We'll check if the user wants to continue after preprocessing.
    # If not, we'll call sys.exit(0) as before. If yes, we'll just return.
    def maybe_exit():
        """Exit only if --continue was NOT passed."""
        if not getattr(args, "continue", False):
            sys.exit(0)
        else:
            logger.info("[PREPROCESSING] --continue specified; not exiting => returning to main.")

    # --------------------------------------------------------------
    # IF --egap_gff NOT PASSED => skip coverage-based pipeline,
    # just run unique_gene_id.py on -a, write "preprocessed_best.gtf".
    # --------------------------------------------------------------
    if not getattr(args, "egap_gff", False):
        logger.info("[PREPROCESSING] => Skipping coverage-based pipeline (no --egap_gff).")
        if not args.a or not os.path.isfile(args.a):
            logger.error("[PREPROCESSING] requires -a <some.gtf> when skipping coverage pipeline.")
            sys.exit(1)

        # Step 2: Download + chmod unique_gene_id.py
        cmd_wget_uniq = (
            "wget -O unique_gene_id.py "
            "https://raw.githubusercontent.com/cfarkas/amyg/refs/heads/main/third_parties/unique_gene_id.py"
        )
        run_pipeline_command(cmd_wget_uniq, args.use_conda, args.use_docker, outdir)
        run_pipeline_command("chmod 755 unique_gene_id.py", args.use_conda, args.use_docker, outdir)

        # Run unique_gene_id.py on user-provided -a
        a_basename = os.path.basename(args.a)
        a_base_noext = os.path.splitext(args.a)[0]
        out_uniq = a_base_noext + ".unique_gene_id.gtf"

        cmd_uniq = f"python unique_gene_id.py {a_basename}"
        logger.info(f"[PREPROCESSING] => {cmd_uniq}")
        run_pipeline_command(cmd_uniq, args.use_conda, args.use_docker, outdir)

        if not os.path.isfile(out_uniq):
            logger.error(f"[PREPROCESSING] missing {out_uniq} => cannot proceed.")
            sys.exit(1)

        # Copy result to "preprocessed_best.gtf"
        final_tuned = os.path.join(outdir, "preprocessed_best.gtf")
        shutil.copy(out_uniq, final_tuned)

        logger.info(f"[PREPROCESSING] => Overriding -a => {final_tuned}")
        args.a = final_tuned

        maybe_exit()  # exit or return, depending on --continue
        return  # in case maybe_exit() does not exit (if --continue was set)

    # --------------------------------------------------------------------
    # ELSE => user passed --egap_gff => run the ENTIRE coverage-based code
    # (All of your original logic remains below, with the -G fix)
    # --------------------------------------------------------------------

    if not args.bam:
        logger.error("[PREPROCESSING] requires --bam <reads.bam>")
        sys.exit(1)
    if not args.egap_gff:
        logger.error("[PREPROCESSING] requires --egap_gff <some.gff>")
        sys.exit(1)

    # Step A: detect coverage
    coverage = detect_coverage(args.bam, outdir, args.use_conda, args.use_docker)
    if coverage is None:
        coverage = 10.0  # fallback if coverage is undetectable

    # Example coverage factors (including 0.125 now)
    factors = [0.125, 0.25, 0.5, 1, 2, 4]
    param_candidates = []
    for fct in factors:
        cval = max(1.0, coverage * fct)
        param_candidates.append(cval)

    logger.info(f"[PREPROCESSING] coverage param sets => {param_candidates}")
    tuned_gtfs = []
    bam_abs = os.path.abspath(args.bam)

    # Step 1: run stringtie with each coverage factor
    for cval in param_candidates:
        lbl = f"cov{cval:.2f}"
        out_gtf = os.path.join(outdir, f"tune_{lbl}.gtf")
        cmd_st = (
            f"stringtie {bam_abs} "
            f"-p {args.threads} "
            f"-c {cval:.2f} "
            f"-s 8 "
            f"-f 0.05 "
            f"-j 3 "
            f"-m 300 "
            f"-o {os.path.basename(out_gtf)}"
        )
        # If user passed --egap_gff, we add -G <annotation> for reference-guided assembly
        if args.egap_gff:
            cmd_st += f" -G {os.path.basename(args.egap_gff)}"

        logger.info(f"[PREPROCESSING] => {cmd_st}")
        run_pipeline_command(cmd_st, args.use_conda, args.use_docker, outdir)
        if os.path.isfile(out_gtf):
            tuned_gtfs.append(out_gtf)

    # Also consider user-provided GTF as a candidate
    if args.a and os.path.isfile(args.a):
        tuned_gtfs.append(args.a)

    # Step 2: Download and chmod unique_gene_id.py
    cmd_wget_uniq = (
        "wget -O unique_gene_id.py "
        "https://raw.githubusercontent.com/cfarkas/amyg/refs/heads/main/third_parties/unique_gene_id.py"
    )
    run_pipeline_command(cmd_wget_uniq, args.use_conda, args.use_docker, outdir)
    run_pipeline_command("chmod 755 unique_gene_id.py", args.use_conda, args.use_docker, outdir)

    best_f1 = -1.0
    best_gtf = None

    # Step 3: For each candidate GTF => run unique_gene_id => run gffcompare => parse F1
    for cand_gtf in tuned_gtfs:
        lbl = os.path.basename(cand_gtf).replace(".gtf", "")
        cand_gtf_base = os.path.splitext(cand_gtf)[0]
        out_uniq = cand_gtf_base + ".unique_gene_id.gtf"

        cmd_uniq = f"python unique_gene_id.py {os.path.basename(cand_gtf)}"
        logger.info(f"[PREPROCESSING] => {cmd_uniq}")
        run_pipeline_command(cmd_uniq, args.use_conda, args.use_docker, outdir)

        if not os.path.isfile(out_uniq):
            alt = glob.glob(cand_gtf_base + ".unique_gene_id.gtf")
            if alt:
                out_uniq = alt[0]
            else:
                logger.warning(f"[PREPROCESSING] missing {out_uniq} => skip")
                continue

        cmp_prefix = f"cmp_{lbl}"
        cmd_cmp = (
            f"gffcompare -r {os.path.basename(args.egap_gff)} "
            f"-o {cmp_prefix}.stats "
            f"{os.path.basename(out_uniq)}"
        )
        logger.info(f"[PREPROCESSING] => {cmd_cmp}")
        run_pipeline_command(cmd_cmp, args.use_conda, args.use_docker, outdir)

        # Fix: rename the stats file if needed
        original_stats = os.path.join(
            outdir, f"{cmp_prefix}.{os.path.basename(out_uniq)}.stats"
        )
        renamed_stats = os.path.join(outdir, f"{cmp_prefix}.stats")
        if os.path.exists(original_stats):
            os.rename(original_stats, renamed_stats)

        stats_file = renamed_stats
        f1_score = 0.0
        if os.path.isfile(stats_file):
            f1_score = compute_f1_from_stats(stats_file)
        else:
            logger.warning(f"compute_f1_from_stats => missing {stats_file}")

        logger.info(f"[PREPROCESSING] => {lbl} => F1={f1_score:.2f}")

        if f1_score > best_f1:
            best_f1 = f1_score
            best_gtf = out_uniq

    if not best_gtf:
        logger.error("[PREPROCESSING] Could not find best GTF => exiting.")
        sys.exit(1)

    # Step 4: Copy the best GTF to "preprocessed_best.gtf" if different
    final_tuned = os.path.join(outdir, "preprocessed_best.gtf")
    if not os.path.exists(final_tuned):
        shutil.copy(best_gtf, final_tuned)
    elif not os.path.samefile(best_gtf, final_tuned):
        shutil.copy(best_gtf, final_tuned)

    # Step 5: Optionally run merge_stringtie_names.py if user passed --egap_gff
    if args.egap_gff:
        logger.info("=== Preprocessing: Downloading merge_stringtie_names.py ===")
        run_pipeline_command(
            "wget https://raw.githubusercontent.com/cfarkas/amyg/refs/heads/main/scripts/merge_stringtie_names.py",
            args.use_conda,
            args.use_docker,
            outdir
        )
        run_pipeline_command("chmod 755 merge_stringtie_names.py", args.use_conda, args.use_docker, outdir)
        out_gtf = "transcripts_named.gtf"
        cmd_merge = (
            f"python merge_stringtie_names.py "
            f"--stringtie_gtf {final_tuned} "
            f"--egap_gff {args.egap_gff} "
            f"--output_gtf {out_gtf}"
        )
        logger.info(f"Running: {cmd_merge}")
        run_pipeline_command(cmd_merge, args.use_conda, args.use_docker, outdir)
        logger.info(f"Transcripts named => {out_gtf}")
        final_tuned = os.path.join(outdir, out_gtf)

    logger.info(f"[PREPROCESSING] best => {best_gtf}, F1={best_f1:.2f}")
    logger.info(f"[PREPROCESSING] => Overriding -a => {final_tuned}")
    args.a = final_tuned

    maybe_exit()  # exit or return, depending on --continue
    # If maybe_exit doesn't exit, we return to main => continue pipeline
    return



###############################################################################
# SINGLE-CELL LOGIC
###############################################################################
def run_stringtie_single_cell(bamfile, out_gtf, ref_gtf, threads):
    cmd=[
        "stringtie",
        "-p", str(threads),
        "-G", os.path.abspath(ref_gtf),
        "-c","2","-s","8","-f","0.05","-j","3","-m","300","-v",
        "-o", out_gtf,
        os.path.abspath(bamfile)
    ]
    logger.info(f"[SINGLE_CELL] => {' '.join(cmd)}")
    try:
        subprocess.run(cmd, check=True)
        return True
    except subprocess.CalledProcessError:
        logger.warning(f"[SINGLE_CELL] StringTie failed => skip {bamfile}")
        return False

def parse_attrs(a_str):
    d={}
    for chunk in a_str.split(';'):
        chunk=chunk.strip()
        if not chunk: 
            continue
        parts=chunk.split(' ',1)
        if len(parts)<2:
            continue
        k,v=parts
        v=v.strip().strip('"')
        d[k]=v
    return d

def interval_overlap_len(ranges, st, en):
    ov=0
    for (s,e) in ranges:
        if e<st:
            continue
        if s>en:
            break
        ovs=max(s,st)
        ove=min(e,en)
        if ove>=ovs:
            ov+=(ove-ovs+1)
    return ov

def merge_intervals(ivals):
    ivals.sort(key=lambda x:x[0])
    merged=[]
    cs,ce=ivals[0]
    for (s,e) in ivals[1:]:
        if s<=ce+1:
            ce=max(ce,e)
        else:
            merged.append((cs,ce))
            cs,ce=s,e
    merged.append((cs,ce))
    return merged

def load_ref_intervals_for_sc(ref_gtf):
    ref_gtf_abs=os.path.abspath(ref_gtf)
    known={}
    with open(ref_gtf_abs) as rf:
        for line in rf:
            if line.startswith('#'):
                continue
            f=line.strip().split('\t')
            if len(f)<5: 
                continue
            c=f[0]
            st=int(f[3])
            en=int(f[4])
            if c not in known:
                known[c]=[]
            known[c].append((st,en))
    for c in known:
        known[c]=merge_intervals(known[c])
    return known

def parse_and_filter_strg(in_gtf, ref_intervals, overlap_thr, library_prefix, merged_lines, exon_info):
    trans_keep={}
    with open(in_gtf) as fin:
        for line in fin:
            if line.startswith('#') or not line.strip():
                continue
            f=line.rstrip().split('\t')
            if len(f)<9:
                continue
            chrom,start,end,feat=f[0],int(f[3]),int(f[4]),f[2]
            attr=parse_attrs(f[8])
            gene_id=attr.get("gene_id","")
            if not gene_id.startswith("STRG."):
                continue
            tid=attr.get("transcript_id","")
            if feat=="transcript":
                length=end-start+1
                if length<=0:
                    continue
                overlap=0
                if chrom in ref_intervals:
                    overlap=interval_overlap_len(ref_intervals[chrom],start,end)
                frac=overlap/length
                keep=(frac<=overlap_thr)
                trans_keep[tid]=keep
                exon_info[tid]={"exons":0,"library":library_prefix}
            elif feat=="exon":
                if tid in trans_keep and trans_keep[tid]:
                    exon_info[tid]["exons"]+=1
    keep_count=0
    all_lines=[]
    with open(in_gtf) as fin:
        for l in fin:
            if l.strip() and not l.startswith('#'):
                all_lines.append(l.rstrip('\n'))
    for line in all_lines:
        f=line.split('\t')
        if len(f)<9:
            continue
        attr=parse_attrs(f[8])
        gene_id=attr.get("gene_id","")
        tid=attr.get("transcript_id","")
        if gene_id.startswith("STRG.") and tid in trans_keep:
            if trans_keep[tid]:
                merged_lines.append(line+'\n')
                keep_count+=1
    return keep_count

def parse_final_annot_sc(final_gtf, exon_info, outdir):
    transcript_map={}
    with open(final_gtf) as fin:
        for line in fin:
            if line.startswith('#') or not line.strip():
                continue
            f=line.rstrip().split('\t')
            if len(f)<9:
                continue
            attr=parse_attrs(f[8])
            tid=attr.get("transcript_id","")
            if tid not in transcript_map:
                transcript_map[tid]={"blastx_hit":"no_hits","gawn_name":"no_hits"}
            if "blastx_hit" in attr:
                transcript_map[tid]["blastx_hit"]=attr["blastx_hit"]
            if "gawn_name" in attr:
                transcript_map[tid]["gawn_name"]=attr["gawn_name"]
    csv_path=os.path.join(outdir,"final_report.csv")
    with open(csv_path,"w") as rep:
        rep.write("transcript_id,library_id,exon_count,blastx_hit,gawn_name\n")
        for tid,data in exon_info.items():
            ex_count=data["exons"]
            lib_id=data["library"]
            b=transcript_map.get(tid,{}).get("blastx_hit","no_hits")
            g=transcript_map.get(tid,{}).get("gawn_name","no_hits")
            rep.write(f"{tid},{lib_id},{ex_count},{b},{g}\n")

    multi_gtf=os.path.join(outdir,"multiexon_only.gtf")
    lines=[]
    with open(final_gtf) as fin:
        for l in fin:
            if l.strip() and not l.startswith('#'):
                lines.append(l.rstrip('\n'))
    multi=[]
    for ln in lines:
        f=ln.split('\t')
        if len(f)<9:
            continue
        attr=parse_attrs(f[8])
        tid=attr.get("transcript_id","")
        if tid in exon_info:
            if exon_info[tid]["exons"]>1:
                multi.append(ln)
    with open(multi_gtf,"w") as mout:
        mout.write("# multi-exon only transcripts\n")
        for line in multi:
            mout.write(line+'\n')

def run_single_cell_mode(args):
    date_str=datetime.datetime.now().strftime("%Y%m%d")
    outdir=os.path.join(args.output, f"amyg_singlecell_{date_str}")
    os.makedirs(outdir, exist_ok=True)

    if not args.input_dir or not args.ref_gtf or not args.ref_fa:
        print("[ERROR] single_cell mode => requires --input_dir, --ref_gtf, --ref_fa", file=sys.stderr)
        sys.exit(1)
    ref_gtf_abs = os.path.abspath(args.ref_gtf)
    ref_fa_abs  = os.path.abspath(args.ref_fa)
    bam_list=glob.glob(os.path.join(os.path.abspath(args.input_dir),"*.bam"))
    if not bam_list:
        print("[ERROR] no .bam found => abort", file=sys.stderr)
        sys.exit(1)

    known={}
    if os.path.isfile(ref_gtf_abs):
        known=load_ref_intervals_for_sc(ref_gtf_abs)

    merged_novel=[]
    exon_info={}
    for bamfile in bam_list:
        prefix=os.path.splitext(os.path.basename(bamfile))[0]
        out_gtf=os.path.join(outdir, f"{prefix}.stringtie.gtf")
        ok=run_stringtie_single_cell(bamfile, out_gtf, ref_gtf_abs, args.threads)
        if not ok:
            continue
        kept=parse_and_filter_strg(out_gtf, known, args.overlap_frac, prefix, merged_novel,exon_info)
        if kept<1:
            print(f"[WARN] no novel => skip {prefix}")

    if not merged_novel:
        print("[WARN] no novel => skip final pipeline.")
        sys.exit(0)

    all_novel=os.path.join(outdir,"all_novel.gtf")
    with open(all_novel,"w") as an:
        an.write("# all novel transcripts\n")
        for l in merged_novel:
            an.write(l)
    print("[INFO] final pass => normal pipeline => all_novel.gtf")
    final_merged_dir=os.path.join(outdir,"merged_final")
    os.makedirs(final_merged_dir, exist_ok=True)

    cmd_amyg=[
        "python",
        os.path.abspath(sys.argv[0]),
        "-a", os.path.abspath(all_novel),
        "-g", ref_fa_abs,
        "--threads", str(args.threads),
        "-o", os.path.abspath(final_merged_dir),
        "--force",
        "--use_docker"
    ]
    try:
        subprocess.run(cmd_amyg, check=True)
    except subprocess.CalledProcessError:
        print("[WARN] final pipeline pass failed, continuing anyway")

    final_annot=os.path.join(final_merged_dir,"final_results","final_annotated.gtf")
    if not os.path.isfile(final_annot):
        print("[WARN] no final_annotated.gtf => can't produce final report")
        sys.exit(0)
    print("[INFO] generating final_report.csv, multiexon_only.gtf => final_results")
    parse_final_annot_sc(final_annot,exon_info, os.path.join(final_merged_dir,"final_results"))
    print("[INFO] single_cell pipeline => done")
    sys.exit(0)

###############################################################################
# MAIN
###############################################################################
def main():
    parser = argparse.ArgumentParser(
        description="annotation pipeline with optional single_cell or coverage-based preprocessing."
    )
    parser.add_argument("--install", choices=["conda", "docker"], help="Install environment and exit.")
    parser.add_argument("--use_conda", action="store_true", help="Use conda env.")
    parser.add_argument("--use_docker", action="store_true", help="Use docker image.")
    parser.add_argument("--threads", type=int, default=10)
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--purge_all_envs", action="store_true")
    parser.add_argument("--dups", action="store_true")
    parser.add_argument("--chunk_size", type=int, default=20000)
    parser.add_argument("-o", "--output", help="Output directory")
    parser.add_argument("-a", help="GTF from StringTie or param tuning")
    parser.add_argument("-g", help="Reference genome (FASTA)")
    parser.add_argument("--egap_gff", help="EGAP GFF for merging or reference checks")
    parser.add_argument("--single_cell", action="store_true",
                        help="If set => multi-bam single-cell logic => exit after done.")
    parser.add_argument("--input_dir", help="Directory with .bam for single_cell")
    parser.add_argument("--ref_gtf", help="Known ref GTF for single_cell mode")
    parser.add_argument("--ref_fa", help="Reference genome for single_cell final pass")
    parser.add_argument("--overlap_frac", type=float, default=0.05,
                        help="Overlap fraction for single_cell mode.")
    parser.add_argument("--preprocessing", action="store_true",
                        help="If set => coverage-based param sets => pick best => override -a.")
    parser.add_argument("--bam", help="BAM used for coverage detection in preprocessing")

    # NEW FLAG
    parser.add_argument("--continue", action="store_true",
                        help="If set, continue the pipeline after preprocessing instead of exiting.")

    args = parser.parse_args()

    if args.purge_all_envs:
        purge_all_envs()
    if args.install == "conda":
        install_conda_env()
    elif args.install == "docker":
        install_docker_image()

    # NEW: unify absolute paths and local copies
    check_inputs(args)

    if args.single_cell:
        run_single_cell_mode(args)

    if args.preprocessing:
        # run_preprocessing_tuning now checks args.continue
        # If not set, it calls sys.exit(0). If set, it returns normally
        run_preprocessing_tuning(args)

    # We only reach here if:
    #   1) --preprocessing wasn't used, OR
    #   2) --preprocessing was used with --continue, so we returned
    log_green_info("Starting normal pipeline script...")

    if not args.a or not args.g or not args.output:
        logger.error("For normal pipeline => must provide -a, -g, -o.")
        sys.exit(1)

    use_conda = args.use_conda
    use_docker = args.use_docker
    threads = args.threads
    force = args.force
    dups = args.dups
    chunk_size = args.chunk_size
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
            gc.write('SWISSPROT_DB="TO_BE_REPLACED"\n')
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
            if "SWISSPROT_DB" in line and "TO_BE_REPLACED" in line:
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

    if args.dups:
        logger.info("::: Step 9B => duplicates => amyg_syntenyblast.py :::")
        run_pipeline_command(
            "curl -O https://raw.githubusercontent.com/cfarkas/amyg/refs/heads/main/scripts/amyg_syntenyblast.py",
            use_conda, use_docker, output_dir
        )
        run_pipeline_command("chmod +x amyg_syntenyblast.py", use_conda, use_docker, output_dir)
        local_output_dir = "."
        if use_docker:
            cmd_synteny = (
                f"PYTHONUNBUFFERED=1 python amyg_syntenyblast.py "
                f"--fasta {os.path.basename(g_filename)} "
                f"--output_dir {local_output_dir} "
                f"--chunk_size {chunk_size} "
                f"--threads {threads}"
            )
        else:
            cmd_synteny = (
                f"python amyg_syntenyblast.py "
                f"--fasta {os.path.basename(g_filename)} "
                f"--output_dir {local_output_dir} "
                f"--chunk_size {chunk_size} "
                f"--threads {threads}"
            )
        run_pipeline_command(cmd_synteny, use_conda, use_docker, output_dir)
        logger.info("::: duplication analysis => done :::")

        logger.info("::: Step 9C => amyg_annotatedups.py for GTF duplication annotation :::")
        run_pipeline_command(
            "curl -O https://raw.githubusercontent.com/cfarkas/amyg/refs/heads/main/scripts/amyg_annotatedups.py",
            use_conda, use_docker, output_dir
        )
        run_pipeline_command("chmod +x amyg_annotatedups.py", use_conda, use_docker, output_dir)

        if use_docker:
            final_annotated_gtf_path = "/data/final_results/final_annotated.gtf"
            synteny_csv_path = "/data/synteny_blocks.csv"
            final_annot_dups_path = "/data/final_results/final_annotated_dups.gtf"
            dup_annot_log = "/data/dup_annot.log"
        else:
            final_annotated_gtf_path = os.path.join(output_dir, "final_results", "final_annotated.gtf")
            synteny_csv_path = os.path.join(".", "synteny_blocks.csv")
            final_annot_dups_path = os.path.join(output_dir, "final_results", "final_annotated_dups.gtf")
            dup_annot_log = os.path.join(".", "dup_annot.log")

        annotate_dups_cmd = (
            f"python amyg_annotatedups.py "
            f"{final_annotated_gtf_path} "
            f"{synteny_csv_path} "
            f"{final_annot_dups_path} "
            f"{dup_annot_log}"
        )
        run_pipeline_command(annotate_dups_cmd, use_conda, use_docker, output_dir)
        logger.info("::: duplication annotation => done :::")

        synteny_csv_on_host = os.path.join(output_dir, "synteny_blocks.csv")
        if os.path.isfile(synteny_csv_on_host):
            shutil.move(synteny_csv_on_host, os.path.join(output_dir, "final_results"))
        pdfs = glob.glob(os.path.join(output_dir, "*.pdf"))
        for pdf_file in pdfs:
            shutil.move(pdf_file, os.path.join(output_dir, "final_results"))

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
