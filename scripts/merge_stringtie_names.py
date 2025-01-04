#!/usr/bin/env python3

"""
merge_stringtie_names_2pass.py

Implements a pipeline to:
 1) Filter both GTFs by 'transcript_id' 
 2) Check for 'gene_id' presence
 3) Fix them with gffread
 4) Run gffcompare with -p prefix + wildcard detection
 5) Find the *annotated.gtf
 6) TWO-PASS approach to unify gene_id across transcripts/exons:
    - Pass 1: parse transcript lines => transcript_id -> final_gene_id 
        * If gene_name != '.' => final_gene_id = gene_name
        * Else final_gene_id = that transcript line's gene_id
    - Pass 2: rewrite all lines (transcript + exon), 
        * Overwrite gene_id with final_gene_id from the dictionary

This ensures exons match the transcript-level gene_id from the final annotation.

Usage:
  python merge_stringtie_names_2pass.py \
    --stringtie_gtf /path/to/stringtie.gtf \
    --egap_gtf /path/to/reference.gtf \
    --prefix gffcmp_out \
    --output_annotated unified_annotated.gtf

Dependencies:
  - Python >=3
  - 'gffread', 'gffcompare' in PATH
"""

import argparse
import os
import sys
import glob
import subprocess

def parse_args():
    p = argparse.ArgumentParser(description="Two-pass unify gene_id in annotated GTF from gffcompare.")
    p.add_argument("--stringtie_gtf", required=True,
                   help="Path to the StringTie GTF")
    p.add_argument("--egap_gtf", required=True,
                   help="Path to the NCBI EGAP GTF (reference).")
    p.add_argument("--prefix", default="gffcmp_out",
                   help="Prefix for gffcompare (-o prefix -p prefix). Default=gffcmp_out")
    p.add_argument("--output_annotated", default="annotated_with_renamed.gtf",
                   help="Final annotated GTF after 2-pass gene_id unify. Default=annotated_with_renamed.gtf")
    return p.parse_args()

def run_cmd(cmd_list, msg=""):
    print(f"[INFO] {msg} Running:\n   {' '.join(cmd_list)}\n")
    cp = subprocess.run(cmd_list)
    if cp.returncode != 0:
        print(f"[ERROR] Command failed (exit code {cp.returncode}): {' '.join(cmd_list)}", file=sys.stderr)
        sys.exit(cp.returncode)

def filter_gtf(in_gtf, out_gtf, label):
    """
    Keep only lines containing 'transcript_id'.
    Print how many lines are kept/skipped, plus examples.
    """
    print(f"[INFO] Filtering lines w/o transcript_id in {in_gtf} => {out_gtf} ({label})")
    kept = 0
    skipped = 0
    first_skips = []
    with open(in_gtf) as fin, open(out_gtf, "w") as fout:
        for line in fin:
            if line.startswith("#"):
                fout.write(line)
                continue
            if "transcript_id" in line:
                fout.write(line)
                kept += 1
            else:
                skipped += 1
                if len(first_skips) < 5:
                    first_skips.append(line.rstrip("\n"))
    print(f"[INFO] {label}: Kept {kept} lines; Skipped {skipped} lines.")
    if first_skips:
        print(f"[INFO] Example skipped lines ({label}):")
        for s in first_skips:
            print("   ", s)
    print("")
    return out_gtf

def check_for_missing_gene_id(in_gtf, label):
    """
    Check how many lines are missing gene_id or transcript_id,
    up to 5 examples.
    """
    missing_count = 0
    examples = []
    with open(in_gtf) as f:
        for idx, line in enumerate(f, start=1):
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            attr = parts[8]
            if ("gene_id" not in attr) or ("transcript_id" not in attr):
                missing_count += 1
                if len(examples)<5:
                    examples.append((idx, line.strip()))
    if missing_count>0:
        print(f"[DEBUG] {label}: Found {missing_count} lines missing gene_id or transcript_id. Examples:")
        for i, txt in examples:
            print(f"   line {i}: {txt}")
    else:
        print(f"[DEBUG] {label}: All lines appear to have gene_id + transcript_id (or comment).")

def gffread_fix(in_gtf, out_gtf, label):
    """
    gffread in_gtf -T -o out_gtf
    """
    cmd = ["gffread", in_gtf, "-T", "-o", out_gtf]
    run_cmd(cmd, f"[{label}] gffread fix")
    if not os.path.isfile(out_gtf):
        print(f"[ERROR] gffread did not produce {out_gtf} for {label}", file=sys.stderr)
        sys.exit(1)
    return out_gtf

def run_gffcompare(sty_fixed, ref_fixed, prefix):
    """
    gffcompare -r ref_fixed -o prefix -p prefix sty_fixed
    => produce prefix*.annotated.gtf, prefix*.tmap, etc.
    """
    cmd = ["gffcompare", "-r", ref_fixed, "-o", prefix, "-p", prefix, sty_fixed]
    run_cmd(cmd, "[gffcompare]")

def parse_attr(attr_str):
    """
    Convert GTF attribute string => dict: key->val
    e.g. gene_id "STRG.12"; gene_name "CLC2DL3";
    """
    d = {}
    for item in attr_str.strip().strip(";").split(";"):
        item = item.strip()
        if not item:
            continue
        parts = item.split(" ",1)
        if len(parts)==2:
            k = parts[0]
            v = parts[1].strip().strip('"')
            d[k] = v
    return d

def build_attr(d):
    items = []
    for k,v in d.items():
        items.append(f'{k} "{v}"')
    return "; ".join(items) + ";"

def pass1_build_transcript_map(in_gtf):
    """
    Pass 1: read all lines, but only feature=transcript lines
    or lines that have 'transcript_id' with possibly gene_name.
    We store: transcript_id -> final_gene_id, where final_gene_id
    = gene_name if gene_name!='.', else = the line's gene_id.
    Return dict tmap: { transcript_id -> final_gene_id }
    """
    tmap = {}
    with open(in_gtf) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts)<9:
                continue
            feature = parts[2]
            # We specifically want lines describing 'transcript'
            # but you can also do: if 'transcript_id' in parts[8] ...
            if feature == "transcript":
                attr_d = parse_attr(parts[8])
                if "transcript_id" not in attr_d:
                    continue
                tid = attr_d["transcript_id"]
                # pick finalGene
                if "gene_name" in attr_d and attr_d["gene_name"] != ".":
                    finalGene = attr_d["gene_name"]
                else:
                    finalGene = attr_d.get("gene_id","?")
                tmap[tid] = finalGene
    return tmap

def pass2_rewrite_unify(in_gtf, out_gtf, transcript_map):
    """
    Pass 2: rewrite all lines, for any that have transcript_id in transcript_map,
    we unify gene_id = transcript_map[that tid].
    We preserve all other attributes.
    """
    changed = 0
    with open(in_gtf) as fin, open(out_gtf, "w") as fout:
        for line in fin:
            if line.startswith("#"):
                fout.write(line)
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts)<9:
                fout.write(line)
                continue

            attr_d = parse_attr(parts[8])
            tid = attr_d.get("transcript_id", None)
            if tid and tid in transcript_map:
                newGene = transcript_map[tid]
                # Overwrite
                if "gene_id" in attr_d:
                    old = attr_d["gene_id"]
                    if old!=newGene:
                        attr_d["gene_id"] = newGene
                        changed += 1
                else:
                    # no gene_id => create it
                    attr_d["gene_id"] = newGene
                    changed += 1

            parts[8] = build_attr(attr_d)
            fout.write("\t".join(parts) + "\n")
    print(f"[INFO] pass2_rewrite_unify: Overwrote gene_id in {changed} lines total.\n")

def main():
    args = parse_args()

    # Decide final annotated
    out_annot = os.path.abspath(args.output_annotated)
    out_dir = os.path.dirname(out_annot)
    if not os.path.isdir(out_dir):
        print(f"[ERROR] Output directory {out_dir} does not exist!")
        sys.exit(1)

    # 1) Filter input GTFs
    sty_filter = os.path.join(out_dir, "stringtie_filter.gtf")
    eg_filter  = os.path.join(out_dir, "egap_filter.gtf")
    filter_gtf(args.stringtie_gtf, sty_filter, "StringTie")
    filter_gtf(args.egap_gtf,      eg_filter,  "EGAP")

    # 2) Check for gene_id
    print("[INFO] Checking gene_id+transcript_id presence...\n")
    check_for_missing_gene_id(sty_filter, "StringTie (filtered)")
    check_for_missing_gene_id(eg_filter,  "EGAP (filtered)")

    # 3) gffread fix => sty_fixed, eg_fixed
    sty_fixed = os.path.join(out_dir, "stringtie_fixed.gtf")
    eg_fixed  = os.path.join(out_dir, "egap_fixed.gtf")
    gffread_fix(sty_filter, sty_fixed, "StringTie")
    gffread_fix(eg_filter,  eg_fixed,  "EGAP")

    # 4) run gffcompare => produce prefix*.annotated.gtf
    print(f"[INFO] Running gffcompare => prefix={args.prefix}\n")
    run_gffcompare(sty_fixed, eg_fixed, args.prefix)

    # 5) find the annotated gtf
    ann_pattern = args.prefix + "*.annotated.gtf"
    ann_matches = glob.glob(ann_pattern)
    if len(ann_matches)==1:
        annotated_gtf = ann_matches[0]
        print(f"[INFO] Found annotated GTF => {annotated_gtf}\n")
    else:
        print(f"[ERROR] Did not uniquely find annotated GTF among => {ann_matches}")
        sys.exit(1)

    # 6) Two-pass unify gene_id
    print("[INFO] Pass 1: building transcript->finalGene map from annotated GTF (transcript lines only).")
    tmap = pass1_build_transcript_map(annotated_gtf)
    print(f"[INFO] Found transcript_id -> finalGene for {len(tmap)} transcripts.\n")

    print("[INFO] Pass 2: rewriting entire annotated GTF => unify exons with that final gene_id.\n")
    pass2_rewrite_unify(annotated_gtf, out_annot, tmap)

    print("[INFO] Done.\n")
    print(f"Final => {out_annot}")
    print("Intermediate files in =>", out_dir)

if __name__=="__main__":
    main()
