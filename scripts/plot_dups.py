#!/usr/bin/env python3

"""
plot_dups_split_all_duplications.py

Generates TWO plots:

Plot A:
  - Ancient (green), Recent (orange), Other (gray)
  - Overall duplication% in upper-left corner
  - Legend outside (frame off)
  - Called "[Plot A]" in logs

Plot B:
  - Classify *all duplications* (either “ancient” or “recent”) as:
      * "intra-only" (blue)   if contig is found ONLY in self-synteny blocks
      * "inter-only" (red)    if contig is found ONLY in cross-synteny blocks
      * "both"       (black)  if contig is found in both self- & cross-synteny
      * "other"      (gray)   if duplication_type ∉ {ancient,recent} 
                              OR contig not found in synteny at all
    So the bar shows how we partition the entire set of duplicated genes 
    among (intra-only, inter-only, both). 
  - Upper-left annotation:
       Duplication: 24.78%
       Of duplicated:
         X% intra-only,
         Y% inter-only,
         Z% both
  - X+Y+Z = ~100% of all duplicated genes (ancient+recent).

Usage:
  python plot_dups_split_all_duplications.py \
    -a <annotation.tsv> -g <gtf_file> -s <synteny.csv> -o <outdir> [--num_contigs=<int>]

Default: --num_contigs=100 => hide x-axis if # contigs > 100
"""

import sys
import getopt
import os
import re

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

#############################
# Parse CLI
#############################

def parse_arguments(argv):
    annotation_file = None
    gtf_file = None
    synteny_file = None
    outdir = "./"
    contig_threshold = 100

    try:
        opts, args = getopt.getopt(
            argv,
            "a:g:s:o:",
            ["annotation=", "gtf=", "synteny=", "outdir=", "num_contigs="]
        )
    except getopt.GetoptError as e:
        print(e)
        print("Usage: plot_dups_split_all_duplications.py "
              "-a <annotation> -g <gtf> -s <synteny.csv> -o <outdir> [--num_contigs=<int>]")
        sys.exit(2)

    for opt, arg in opts:
        if opt in ("-a", "--annotation"):
            annotation_file = arg
        elif opt in ("-g", "--gtf"):
            gtf_file = arg
        elif opt in ("-s", "--synteny"):
            synteny_file = arg
        elif opt in ("-o", "--outdir"):
            outdir = arg
        elif opt == "--num_contigs":
            contig_threshold = int(arg)

    if not (annotation_file and gtf_file and synteny_file):
        print("ERROR: must specify -a, -g, and -s.")
        sys.exit(2)

    return {
        "annotation_file": annotation_file,
        "gtf_file": gtf_file,
        "synteny_file": synteny_file,
        "outdir": outdir,
        "contig_threshold": contig_threshold
    }


#############################
# GTF parser
#############################

def parse_gtf_to_dataframe(gtf_file):
    rows = []
    with open(gtf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts)<9:
                continue

            seqname, source, feature, start, end, score, strand, frame, attributes = parts
            if feature!="transcript":
                continue

            attr_dict = {}
            for kv in attributes.split(';'):
                kv=kv.strip()
                if kv:
                    kvpair = kv.split(' ',1)
                    if len(kvpair)==2:
                        key=kvpair[0]
                        val=kvpair[1].replace('"','')
                        attr_dict[key] = val
            t_id = attr_dict.get("transcript_id")
            dup_type = attr_dict.get("duplication_type")

            if t_id:
                rows.append({
                    "Name": t_id,
                    "duplication_type": dup_type,
                    "start": int(start),
                    "end": int(end),
                    "strand": strand,
                    "seqname": seqname
                })
    return pd.DataFrame(rows)


#############################
# Synteny parser
#############################

def parse_synteny_blocks(synteny_file):
    df_syn = pd.read_csv(synteny_file, sep=',', dtype=str).fillna('')
    return df_syn


#############################
# Plot A: ancient/recent/other
#############################

def plot_original_stacked_bar(merged_df, outdir, contig_threshold, dup_percent):
    """
    [Plot A]:
      - ancient=green, recent=orange, other=gray
      - duplication% in top-left
      - legend outside frame
    """
    group_data=[]
    for seq, grp in merged_df.groupby("seqname"):
        if seq=='':
            continue
        c_len = grp["end"].max()-grp["start"].min()+1
        a = (grp["duplication_type"]=="ancient").sum()
        r = (grp["duplication_type"]=="recent").sum()
        o = len(grp) - a - r
        group_data.append({
            "seqname": seq,
            "length": c_len,
            "ancient": a,
            "recent": r,
            "other": o
        })

    dfc = pd.DataFrame(group_data)
    if dfc.empty:
        print("[Plot A] no data found => skipping.")
        return

    dfc = dfc.sort_values("length")
    x_pos = np.arange(len(dfc))
    lens = dfc["length"].astype(float).values
    max_len = lens.max() if len(lens)>0 else 1
    scale_factor = 1.0/max_len
    scaled_w = lens*scale_factor

    anc = dfc["ancient"].values
    rec = dfc["recent"].values
    oth = dfc["other"].values

    plt.figure(figsize=(12,6))
    plt.bar(x_pos, anc, width=scaled_w, color="green", label="ancient")
    plt.bar(x_pos, rec, bottom=anc, width=scaled_w, color="orange", label="recent")
    bottom2 = anc + rec
    plt.bar(x_pos, oth, bottom=bottom2, width=scaled_w, color="gray", label="other")

    if len(dfc)>contig_threshold:
        print(f"[Plot A] {len(dfc)} contigs > {contig_threshold}, hiding x-axis.")
        plt.xticks(x_pos, [""]*len(x_pos))
    else:
        plt.xticks(x_pos, dfc["seqname"], rotation=90)

    # annotate
    plt.text(
        0.01, 0.95,
        f"Duplication: {dup_percent:.2f}%",
        transform=plt.gca().transAxes,
        fontsize=12,
        ha='left', va='top',
        bbox=dict(boxstyle="square,pad=0.3", fc="white", ec="black", alpha=0.7)
    )

    plt.xlabel("Contigs (ordered by ascending length)", fontsize=14)
    plt.ylabel("Number of Genes", fontsize=14)
    plt.title("Stacked Bar: Genes by Duplication Type per Contig", fontsize=14)
    plt.legend(frameon=False, bbox_to_anchor=(1.02,1), loc="upper left")

    outpdf = os.path.join(outdir,"contig_stacked_bar.pdf")
    plt.tight_layout()
    plt.savefig(outpdf, bbox_inches='tight')
    plt.close()
    print(f"[Plot A] => saved {outpdf}")

#############################
# Plot B: split all "ancient" or "recent" into (intra-only, inter-only, both)
#############################

def classify_all_dup_intra_inter(merged_df, synteny_df):
    """
    We'll define "all duplicated genes" => duplication_type in {ancient,recent}.
    For the rest => "other".

    Then for each contig, we see if it shows up in synteny blocks with qc==sc (call that "self") 
    or qc!=sc (call that "cross").
    - If contig is in self only => "intra-only"
    - If contig is in cross only => "inter-only"
    - If contig is in both => "both"
    - If gene is not duplicated => "other"
    - If contig not in synteny => also "other".

    We'll store that classification in "dup_intra_inter_cat" for each gene. 
    Then we compute the fraction of duplicated genes that are "intra-only", "inter-only", "both".
    """

    # gather sets of contigs from synteny
    contigs_self = set()
    contigs_cross= set()
    for idx, row in synteny_df.iterrows():
        qc = row["query_contig"]
        sc = row["subject_contig"]
        if qc=='' or sc=='':
            continue
        if qc==sc:
            contigs_self.add(qc)
        else:
            contigs_cross.add(qc)
            contigs_cross.add(sc)

    cat_list = []
    for idx, row in merged_df.iterrows():
        dt = row["duplication_type"]
        sn = row["seqname"]
        if dt not in ("ancient","recent"):
            # not duplicated => "other"
            cat_list.append("other")
            continue
        # dt is in {ancient,recent} => we classify by contig
        if not sn:
            cat_list.append("other")
            continue
        in_self  = (sn in contigs_self)
        in_cross = (sn in contigs_cross)
        if in_self and not in_cross:
            cat_list.append("intra-only")
        elif not in_self and in_cross:
            cat_list.append("inter-only")
        elif in_self and in_cross:
            cat_list.append("both")
        else:
            # contig not found in synteny => "other"
            cat_list.append("other")

    merged_df["dup_intra_inter_cat"] = cat_list

    # compute fraction among duplicated genes
    dup_mask = merged_df["duplication_type"].isin(["ancient","recent"])
    total_dup = dup_mask.sum()
    if total_dup==0:
        return merged_df, 0.0, 0.0, 0.0

    cat_counts = merged_df.loc[dup_mask,"dup_intra_inter_cat"].value_counts()
    n_intra_only = cat_counts.get("intra-only",0)
    n_inter_only = cat_counts.get("inter-only",0)
    n_both       = cat_counts.get("both",0)

    frac_i = 100.0*n_intra_only/total_dup
    frac_r = 100.0*n_inter_only/total_dup
    frac_b = 100.0*n_both/total_dup
    return merged_df, frac_i, frac_r, frac_b


def plot_intra_inter_stacked_bar(merged_df, outdir, contig_threshold, dup_percent,
                                 frac_intra_only, frac_inter_only, frac_both):
    """
    [Plot B]:
      - "intra-only" => blue
      - "inter-only" => red
      - "both"       => black
      - "other"      => gray (non-duplicated or contig not in synteny)
    """
    group_data = []
    all_seqnames = merged_df["seqname"].dropna().unique()
    for seq in all_seqnames:
        sub = merged_df[merged_df["seqname"]==seq]
        c_len = sub["end"].max()-sub["start"].min()+1

        i_count = (sub["dup_intra_inter_cat"]=="intra-only").sum()
        r_count = (sub["dup_intra_inter_cat"]=="inter-only").sum()
        b_count = (sub["dup_intra_inter_cat"]=="both").sum()
        o_count = (sub["dup_intra_inter_cat"]=="other").sum()

        group_data.append({
            "seqname": seq,
            "length": c_len,
            "intra_only": i_count,
            "inter_only": r_count,
            "both": b_count,
            "other": o_count
        })

    dfc = pd.DataFrame(group_data).sort_values("length")
    if dfc.empty:
        print("[Plot B] no data => skipping.")
        return

    x_pos = np.arange(len(dfc))
    lens  = dfc["length"].astype(float).values
    mx    = lens.max() if len(lens)>0 else 1
    sf    = 1.0/mx
    scaled_w = lens*sf

    i_vals = dfc["intra_only"].values
    r_vals = dfc["inter_only"].values
    b_vals = dfc["both"].values
    o_vals = dfc["other"].values

    plt.figure(figsize=(12,6))
    plt.bar(x_pos, i_vals, width=scaled_w, color="blue", label="intra-only")
    bot1 = i_vals
    plt.bar(x_pos, r_vals, bottom=bot1, width=scaled_w, color="red", label="inter-only")
    bot2 = i_vals+r_vals
    plt.bar(x_pos, b_vals, bottom=bot2, width=scaled_w, color="black", label="both")
    bot3 = i_vals+r_vals+b_vals
    plt.bar(x_pos, o_vals, bottom=bot3, width=scaled_w, color="gray", label="other")

    if len(dfc)>contig_threshold:
        print(f"[Plot B] {len(dfc)} contigs > {contig_threshold}, hiding x-axis.")
        plt.xticks(x_pos, [""]*len(x_pos))
    else:
        plt.xticks(x_pos, dfc["seqname"], rotation=90)

    # annotation
    # "Duplication: 24.78%\nOf duplicated:\n  X% intra-only,\n  Y% inter-only,\n  Z% both"
    text_str = (f"Duplication: {dup_percent:.2f}%\n"
                f"Of duplicated:\n"
                f"  {frac_intra_only:.2f}% intra-only,\n"
                f"  {frac_inter_only:.2f}% inter-only,\n"
                f"  {frac_both:.2f}% both")
    plt.text(
        0.01, 0.95,
        text_str,
        transform=plt.gca().transAxes,
        fontsize=12,
        ha='left', va='top',
        bbox=dict(boxstyle="square,pad=0.3", fc="white", ec="black", alpha=0.7)
    )

    plt.xlabel("Contigs (ordered by ascending length)", fontsize=14)
    plt.ylabel("Number of Genes", fontsize=14)
    plt.title("Plot B: Intra vs Inter among All Duplications", fontsize=14)
    plt.legend(frameon=False, bbox_to_anchor=(1.02,1), loc="upper left")

    outpdf = os.path.join(outdir, "contig_stacked_bar_intra_inter_all_dups.pdf")
    plt.tight_layout()
    plt.savefig(outpdf, bbox_inches='tight')
    plt.close()
    print(f"[Plot B] => saved {outpdf}")


#############################
# main
#############################

def main():
    args = parse_arguments(sys.argv[1:])
    annotation_file = args["annotation_file"]
    gtf_file = args["gtf_file"]
    synteny_file = args["synteny_file"]
    outdir = args["outdir"]
    contig_threshold = args["contig_threshold"]

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    # read annotation + gtf => merged
    annot_df = pd.read_csv(annotation_file, sep='\t', dtype=str).fillna('')
    gtf_df   = parse_gtf_to_dataframe(gtf_file)
    merged_df= pd.merge(annot_df, gtf_df, on="Name", how="left")

    # compute overall duplication% (ancient+recent)
    n_all = len(merged_df)
    dup_mask = merged_df["duplication_type"].isin(["ancient","recent"])
    dup_count= dup_mask.sum()
    dup_percent = (dup_count / n_all * 100) if n_all>0 else 0.0
    print(f"[INFO] {n_all} total transcripts loaded.")
    print(f"[INFO] duplication => {dup_count}/{n_all} = {dup_percent:.2f}%")

    # parse synteny
    syn_df = parse_synteny_blocks(synteny_file)

    # PLOT A
    plot_original_stacked_bar(merged_df, outdir, contig_threshold, dup_percent)

    # PLOT B => classify all duplicates as (intra-only, inter-only, both, other)
    merged_df, frac_i, frac_r, frac_b = classify_all_dup_intra_inter(merged_df, syn_df)
    plot_intra_inter_stacked_bar(
        merged_df,
        outdir,
        contig_threshold,
        dup_percent,
        frac_i, frac_r, frac_b
    )

    # write CSV e.g. "dup_intra_inter.csv" with columns [Name, duplication_type, dup_intra_inter_cat]
    outcsv = os.path.join(outdir, "dup_intra_inter.csv")
    merged_df[["Name","duplication_type","dup_intra_inter_cat"]].to_csv(outcsv, index=False)
    print(f"[Plot B] wrote => {outcsv}")
    print("DONE.")


if __name__=="__main__":
    main()
