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

Additionally:
 7) Add missing genes from EGAP only if their bounding-box overlap
    with existing genes is <= 50%.

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


########################################
# Overlap Checking and Adding New Genes #
########################################

def _compute_gene_bounding_boxes(gtf_path):
    """
    Compute bounding-box coords for each gene in `gtf_path`.
    Return a dict: gene_id -> list of (chr, strand, start, end)
      We store possibly multiple intervals if the GTF has multiple gene lines 
      or complex structure, but often it's just one bounding box per gene.
      However, to keep it simple, we'll unify them into a single min->max
      for each gene_id.

    Note: If the file has no 'gene' feature lines, we fallback by deriving
    bounding boxes from all lines that share the same gene_id.
    """
    from collections import defaultdict

    gene_boxes = {}   # gene_id -> (chr, strand, min_start, max_end)
    stored_chr_strand = {}

    with open(gtf_path, "r") as fin:
        for line in fin:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue

            chrom, source, feature, start, end, score, strand, frame, attr_str = parts
            start = int(start)
            end = int(end)

            attr_d = parse_attr(attr_str)
            gid = attr_d.get("gene_id", None)
            if not gid:
                continue

            if gid not in gene_boxes:
                gene_boxes[gid] = (chrom, strand, start, end)
                stored_chr_strand[gid] = (chrom, strand)
            else:
                # update bounding box
                (c_chrom, c_strand, c_start, c_end) = gene_boxes[gid]
                # If there's a mismatch in chr/strand, it might be an annotation quirk,
                # but typically should not happen. We'll assume they're the same.
                # Just unify min_start, max_end
                gene_boxes[gid] = (
                    c_chrom,
                    c_strand,
                    min(c_start, start),
                    max(c_end, end),
                )
    return gene_boxes


def _compute_overlap_fraction(boxA, boxB):
    """
    Given two bounding boxes (chrA, strandA, startA, endA) and
    (chrB, strandB, startB, endB), compute how much they overlap.
    Return fraction_of_overlap relative to the smaller region.

    If they are on different chromosomes or different strands,
    overlap is 0.

    Overlap fraction = overlap_len / min(lenA, lenB).
    """
    chrA, strandA, startA, endA = boxA
    chrB, strandB, startB, endB = boxB

    if (chrA != chrB) or (strandA != strandB):
        return 0.0

    overlap_start = max(startA, startB)
    overlap_end = min(endA, endB)
    if overlap_end < overlap_start:
        return 0.0  # no overlap
    overlap_len = overlap_end - overlap_start + 1

    lenA = endA - startA + 1
    lenB = endB - startB + 1
    smaller_len = min(lenA, lenB)

    return overlap_len / float(smaller_len)


def add_missing_genes_from_egap(final_gtf, egap_fixed, max_overlap=0.5):
    """
    Checks which gene_ids appear in egap_fixed.gtf but not in final_gtf.
    For any missing gene_id, we check that it does NOT overlap > max_overlap
    fraction with any existing gene. If overlap <= max_overlap, we append.

    Reports how many missing genes were excluded due to overlap > max_overlap
    and how many were appended. Returns the path to an updated GTF file
    that includes the missing lines.
    """
    if not os.path.isfile(final_gtf) or not os.path.isfile(egap_fixed):
        print("[WARN] add_missing_genes_from_egap: one of the input files is missing!")
        return final_gtf  # do nothing

    # Build bounding boxes for existing gene_ids in final_gtf
    existing_boxes = _compute_gene_bounding_boxes(final_gtf)

    # parse gene_ids from final_gtf
    existing_gene_ids = set(existing_boxes.keys())

    # Build bounding boxes for EGAP
    egap_boxes = _compute_gene_bounding_boxes(egap_fixed)

    # parse all lines from egap_fixed, store them by gene_id
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
                egap_map.setdefault(gid, []).append(line)

    # find missing gene_ids
    missing_gene_ids = set(egap_map.keys()) - existing_gene_ids
    if not missing_gene_ids:
        print("[INFO] No missing gene_ids found from EGAP => nothing appended.")
        return final_gtf

    excluded_for_no_box = 0
    excluded_for_overlap = 0
    approved_missing = []

    # Filter out any that overlap > max_overlap fraction with existing boxes
    for mgid in missing_gene_ids:
        if mgid not in egap_boxes:
            # no bounding box found => skip
            excluded_for_no_box += 1
            continue

        box_new = egap_boxes[mgid]
        # check overlap with all existing genes
        overlaps_too_much = False
        for ex_gid, ex_box in existing_boxes.items():
            frac = _compute_overlap_fraction(box_new, ex_box)
            if frac > max_overlap:
                overlaps_too_much = True
                break

        if overlaps_too_much:
            excluded_for_overlap += 1
        else:
            approved_missing.append(mgid)

    if not approved_missing:
        print("[INFO] All missing gene_ids overlapped more than "
              f"{max_overlap*100:.0f}% => none appended.")
        print(f"[INFO] Total missing = {len(missing_gene_ids)}, "
              f"excluded_for_no_box = {excluded_for_no_box}, "
              f"excluded_for_overlap = {excluded_for_overlap}")
        return final_gtf

    # Actually append the approved_missing genes
    out_temp = final_gtf + ".with_egap_missing"
    with open(final_gtf, "r") as oldf, open(out_temp, "w") as newf:
        # copy old lines
        for line in oldf:
            newf.write(line)
        # append missing lines for approved genes
        for mgid in sorted(approved_missing):
            for ln in egap_map[mgid]:
                newf.write(ln)

    os.replace(out_temp, final_gtf)

    print(f"[INFO] add_missing_genes_from_egap: appended lines for {len(approved_missing)} new gene_ids "
          f"(out of {len(missing_gene_ids)} missing).")
    print(f"[INFO]   - excluded_for_no_box = {excluded_for_no_box}")
    print(f"[INFO]   - excluded_for_overlap (> {max_overlap*100:.0f}%) = {excluded_for_overlap}")
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

    # 9) Finally, add missing genes from the EGAP file,
    #    only if overlap <= 50% (set by max_overlap=0.5).
    add_missing_genes_from_egap(out_annot, eg_fixed, max_overlap=0.5)

    print("[INFO] Final attribute reordering + missing gene addition done. See final file =>", out_annot)


if __name__ == "__main__":
    main()
