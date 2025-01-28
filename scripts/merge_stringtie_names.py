#!/usr/bin/env python3

"""
merge_stringtie_names.py

Implements a pipeline to:
 1) Filter both GTF/GFF by 'transcript_id' 
 2) Check for 'gene_id' presence
 3) Fix them with gffread (GFF->GTF conversion if needed)
 4) Run gffcompare with -p prefix + wildcard detection
 5) Find the *annotated.gtf
 6) TWO-PASS approach to unify gene_id across transcripts/exons:
    - Pass 1: parse *transcript lines* => transcript_id -> final_gene_id 
        * If gene_name != '.' => final_gene_id = gene_name
        * Else if cmp_ref_gene != '.' => final_gene_id = cmp_ref_gene
        * Else final_gene_id = that line's gene_id
    - Pass 2: rewrite *every line* (transcript + exon + other features),
        * If it has a transcript_id in our map, unify (overwrite) gene_id with final_gene_id
        * Otherwise leave it as-is

Finally, we reorder the attributes in the final GTF so that
**gene_id** appears first, **transcript_id** second, then everything else.

Usage:
  python merge_stringtie_names.py \
    --stringtie_gtf /path/to/stringtie.gtf \
    --egap_gff /path/to/reference.gff \
    --prefix gffcmp_out \
    --output_gtf transcripts_named.gtf

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
    p = argparse.ArgumentParser(description="Two-pass unify gene_id in annotated GTF from gffcompare, then reorder attributes.")
    p.add_argument("--stringtie_gtf", required=True,
                   help="Path to the StringTie GTF (e.g. output from StringTie).")
    p.add_argument("--egap_gff", required=True,
                   help="Path to the NCBI EGAP GFF (Gnomon-based).")
    p.add_argument("--prefix", default="gffcmp_out",
                   help="Prefix for gffcompare (-o prefix -p prefix). Default=gffcmp_out")
    p.add_argument("--output_annotated", default="annotated_and_renamed.gtf",
                   help="Final annotated GTF after 2-pass gene_id unify. Default=annotated_and_renamed.gtf")

    # (Optional) if you want a secondary alias:
    p.add_argument("--output_gtf",
                   help="Alias for --output_annotated. Use e.g. --output_gtf transcripts_named.gtf")

    args = p.parse_args()

    # If user passes --output_gtf, override args.output_annotated:
    if args.output_gtf:
        args.output_annotated = args.output_gtf

    return args


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
                if len(examples) < 5:
                    examples.append((idx, line.strip()))
    if missing_count > 0:
        print(f"[DEBUG] {label}: Found {missing_count} lines missing gene_id or transcript_id. Examples:")
        for i, txt in examples:
            print(f"   line {i}: {txt}")
    else:
        print(f"[DEBUG] {label}: All lines appear to have gene_id + transcript_id (or comment).")


def gffread_fix(in_gtf, out_gtf, label):
    """
    gffread in_gtf -T -F -o out_gtf
    Using -F if Gnomon-based GFF might be incomplete.
    """
    cmd = ["gffread", in_gtf, "-T", "-F", "-o", out_gtf]
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
    attr_str = attr_str.strip().strip(";")
    for item in attr_str.split(";"):
        item = item.strip()
        if not item:
            continue
        parts = item.split(" ", 1)
        if len(parts) == 2:
            k = parts[0]
            v = parts[1].strip().strip('"')
            d[k] = v
    return d


def build_attr(d):
    """
    Build an attribute string from dict `d`.
    * Reorder so gene_id is first, transcript_id is second, then everything else.
    """
    items = []

    # If gene_id present, put it first
    if "gene_id" in d:
        items.append(f'gene_id "{d["gene_id"]}"')
    # If transcript_id present, put it second
    if "transcript_id" in d:
        items.append(f'transcript_id "{d["transcript_id"]}"')

    # Then the rest
    for k, v in d.items():
        if k not in ["gene_id", "transcript_id"]:
            items.append(f'{k} "{v}"')
    return "; ".join(items) + ";"


def pass1_build_transcript_map(in_gtf):
    """
    Pass 1: read all lines with feature='transcript'.
    transcript_id -> final_gene_id
       If gene_name != '.' => final_gene_id = gene_name
       Else if cmp_ref_gene != '.' => final_gene_id = cmp_ref_gene
       Else => final_gene_id = gene_id
    """
    tmap = {}
    with open(in_gtf) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            feature = parts[2]
            if feature == "transcript":
                attr_d = parse_attr(parts[8])
                if "transcript_id" not in attr_d:
                    continue
                tid = attr_d["transcript_id"]

                # 1) if gene_name present (and not '.'), use it
                if "gene_name" in attr_d and attr_d["gene_name"] != ".":
                    finalGene = attr_d["gene_name"]
                # 2) else if cmp_ref_gene present (and not '.'), use that
                elif "cmp_ref_gene" in attr_d and attr_d["cmp_ref_gene"] != ".":
                    finalGene = attr_d["cmp_ref_gene"]
                # 3) else fallback to gene_id
                else:
                    finalGene = attr_d.get("gene_id", "?")

                tmap[tid] = finalGene
    return tmap


def pass2_rewrite_unify(in_gtf, out_gtf, transcript_map):
    """
    Pass 2: rewrite every line in the annotated GTF,
    unify gene_id for lines that have a transcript_id in transcript_map.

    Then build attributes with build_attr(...) to put gene_id first, then transcript_id, etc.
    """
    changed = 0
    total = 0
    with open(in_gtf) as fin, open(out_gtf, "w") as fout:
        for line in fin:
            if line.startswith("#"):
                fout.write(line)
                continue

            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                fout.write(line)
                continue

            attr_d = parse_attr(parts[8])
            tid = attr_d.get("transcript_id", None)

            if tid and (tid in transcript_map):
                newGene = transcript_map[tid]
                oldGene = attr_d.get("gene_id", None)
                if oldGene != newGene:
                    attr_d["gene_id"] = newGene
                    changed += 1

            parts[8] = build_attr(attr_d)
            fout.write("\t".join(parts) + "\n")
            total += 1

    print(f"[INFO] pass2_rewrite_unify: Overwrote gene_id in {changed} out of {total} lines total.\n")


def remove_prefixes_in_final_gtf(final_gtf_path):
    """
    Additional final step: remove 'gene-' prefix from gene_id/ref_gene_id
    and 'rna-' prefix from transcript_id/cmp_ref, if found.
    Then reorder attributes again (just to be sure).
    """
    tmp_output = final_gtf_path + ".prefix_stripped"
    changed = 0
    total = 0

    with open(final_gtf_path, "r") as fin, open(tmp_output, "w") as fout:
        for line in fin:
            if line.startswith("#"):
                fout.write(line)
                continue

            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                fout.write(line)
                continue

            attr_d = parse_attr(parts[8])

            # remove 'gene-' prefix from gene_id and ref_gene_id
            if "gene_id" in attr_d:
                old_gene = attr_d["gene_id"]
                if old_gene.startswith("gene-"):
                    attr_d["gene_id"] = old_gene.replace("gene-", "", 1)
                    changed += 1

            if "ref_gene_id" in attr_d:
                old_ref_gene = attr_d["ref_gene_id"]
                if old_ref_gene.startswith("gene-"):
                    attr_d["ref_gene_id"] = old_ref_gene.replace("gene-", "", 1)
                    changed += 1

            # remove 'rna-' prefix from transcript_id and cmp_ref
            if "transcript_id" in attr_d:
                old_tx = attr_d["transcript_id"]
                if old_tx.startswith("rna-"):
                    attr_d["transcript_id"] = old_tx.replace("rna-", "", 1)
                    changed += 1

            if "cmp_ref" in attr_d:
                old_cmp = attr_d["cmp_ref"]
                if old_cmp.startswith("rna-"):
                    attr_d["cmp_ref"] = old_cmp.replace("rna-", "", 1)
                    changed += 1

            # re-build attributes w/ gene_id first, transcript_id second
            parts[8] = build_attr(attr_d)
            fout.write("\t".join(parts) + "\n")
            total += 1

    os.replace(tmp_output, final_gtf_path)
    print(f"[INFO] remove_prefixes_in_final_gtf: Stripped prefixes in {changed} out of {total} lines.\n")


def rename_transcript_id_by_cmp_ref(final_gtf_path):
    """
    If a line has 'cmp_ref "<someRef>"' on a transcript, we rename transcript_id => that <someRef>.
    Then apply that rename to all lines with that transcript_id.
    After that, reorder attributes again (just to keep consistent).
    """
    # 1) parse to find transcript feature with transcript_id + cmp_ref => map oldTID -> newTID
    tid_to_cmp = {}
    with open(final_gtf_path, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            if parts[2] == "transcript":
                attr_d = parse_attr(parts[8])
                if "transcript_id" in attr_d and "cmp_ref" in attr_d:
                    old_tid = attr_d["transcript_id"]
                    new_tid = attr_d["cmp_ref"]
                    if new_tid and (new_tid != old_tid):
                        tid_to_cmp[old_tid] = new_tid

    if not tid_to_cmp:
        print("[INFO] rename_transcript_id_by_cmp_ref: no transcripts with cmp_ref found, skipping.\n")
        return

    # 2) rewrite entire GTF => rename transcript_id
    tmp_out = final_gtf_path + ".cmp_ref_renamed"
    changed = 0
    total = 0
    with open(final_gtf_path, "r") as fin, open(tmp_out, "w") as fout:
        for line in fin:
            if line.startswith("#"):
                fout.write(line)
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                fout.write(line)
                continue
            attr_d = parse_attr(parts[8])
            old_tid = attr_d.get("transcript_id", None)
            if old_tid and old_tid in tid_to_cmp:
                new_tid = tid_to_cmp[old_tid]
                attr_d["transcript_id"] = new_tid
                changed += 1

            # reorder
            parts[8] = build_attr(attr_d)
            fout.write("\t".join(parts) + "\n")
            total += 1

    os.replace(tmp_out, final_gtf_path)
    print(f"[INFO] rename_transcript_id_by_cmp_ref: Rewrote transcript_id in {changed} out of {total} lines.\n")


########################
#  ADDED FUNCTION HERE #
########################
def add_missing_genes_from_egap(final_gtf, egap_fixed):
    """
    Checks which gene_ids appear in egap_fixed.gtf but not in final_gtf.
    For any gene_id that is missing, appends all lines from egap_fixed
    with that gene_id. Returns the path to an updated GTF file that includes
    the missing lines, plus prints how many unique gene_ids were appended.
    """
    if not os.path.isfile(final_gtf) or not os.path.isfile(egap_fixed):
        print("[WARN] add_missing_genes_from_egap: one of the input files is missing!")
        return final_gtf  # do nothing

    # parse gene_ids from final_gtf
    existing_gene_ids = set()
    with open(final_gtf, "r") as f_in:
        for line in f_in:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            attr_d = parse_attr(parts[8])
            gid = attr_d.get("gene_id", None)
            if gid:
                existing_gene_ids.add(gid)

    # parse gene_ids from egap_fixed
    # store lines by gene_id
    egap_map = {}
    with open(egap_fixed, "r") as f_in:
        for line in f_in:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            attr_d = parse_attr(parts[8])
            gid = attr_d.get("gene_id", None)
            if gid:
                if gid not in egap_map:
                    egap_map[gid] = []
                egap_map[gid].append(line)

    # find missing gene_ids
    missing_gene_ids = set(egap_map.keys()) - existing_gene_ids

    if not missing_gene_ids:
        print("[INFO] No missing gene_ids found from EGAP => nothing appended.")
        return final_gtf

    # append them to a new file
    out_temp = final_gtf + ".with_egap_missing"
    appended_count = 0
    with open(final_gtf, "r") as oldf, open(out_temp, "w") as newf:
        # first copy old lines
        for line in oldf:
            newf.write(line)
        # now append missing lines
        for mgid in sorted(missing_gene_ids):
            for ln in egap_map[mgid]:
                newf.write(ln)
            appended_count += 1

    os.replace(out_temp, final_gtf)

    print(f"[INFO] add_missing_genes_from_egap: appended lines for {len(missing_gene_ids)} missing gene_ids.")
    return final_gtf


def main():
    args = parse_args()

    # 1) Resolve final annotated output path
    out_annot = os.path.abspath(args.output_annotated)
    out_dir = os.path.dirname(out_annot)
    if not os.path.isdir(out_dir):
        print(f"[ERROR] Output directory {out_dir} does not exist!")
        sys.exit(1)

    # 2) Filter input GTF + GFF to lines that have transcript_id
    sty_filter = os.path.join(out_dir, "stringtie_filter.gtf")
    eg_filter  = os.path.join(out_dir, "egap_filter.gtf")

    print("[INFO] Filtering input files to lines that contain 'transcript_id'...\n")
    filter_gtf(args.stringtie_gtf, sty_filter, "StringTie")
    filter_gtf(args.egap_gff,      eg_filter,  "EGAP GFF")

    # 3) Check for missing gene_id
    print("[INFO] Checking gene_id+transcript_id presence...\n")
    check_for_missing_gene_id(sty_filter, "StringTie (filtered)")
    check_for_missing_gene_id(eg_filter,  "EGAP (filtered)")

    # 4) Use gffread fix => produce sty_fixed, eg_fixed
    sty_fixed = os.path.join(out_dir, "stringtie_fixed.gtf")
    eg_fixed  = os.path.join(out_dir, "egap_fixed.gtf")
    gffread_fix(sty_filter, sty_fixed, "StringTie")
    gffread_fix(eg_filter,  eg_fixed,  "EGAP GFF")

    # 5) run gffcompare => produce prefix*.annotated.gtf
    print(f"[INFO] Running gffcompare => prefix={args.prefix}\n")
    run_gffcompare(sty_fixed, eg_fixed, args.prefix)

    # 6) find the annotated gtf
    ann_pattern = args.prefix + "*.annotated.gtf"
    ann_matches = glob.glob(ann_pattern)
    if len(ann_matches) == 1:
        annotated_gtf = ann_matches[0]
        print(f"[INFO] Found annotated GTF => {annotated_gtf}\n")
    else:
        print(f"[ERROR] Did not uniquely find annotated GTF among => {ann_matches}")
        sys.exit(1)

    # 7) Two-pass unify gene_id
    print("[INFO] Pass 1: building transcript->finalGene map from annotated GTF (transcript lines).")
    tmap = pass1_build_transcript_map(annotated_gtf)
    print(f"[INFO] Found transcript_id -> finalGene for {len(tmap)} transcripts.\n")

    print("[INFO] Pass 2: rewriting entire annotated GTF => unify gene_id...\n")
    pass2_rewrite_unify(annotated_gtf, out_annot, tmap)

    print(f"[INFO] Done pass2 => {out_annot}.\n")

    # 8) Additional final steps:
    #    8a) remove prefixes
    remove_prefixes_in_final_gtf(out_annot)
    #    8b) rename transcript_id by cmp_ref (if found)
    rename_transcript_id_by_cmp_ref(out_annot)

    # 9) Finally, add missing genes from the EGAP file
    #    and see how many new gene_ids are appended:
    add_missing_genes_from_egap(out_annot, eg_fixed)

    print("[INFO] Final attribute reordering done. See final file =>", out_annot)


if __name__ == "__main__":
    main()
