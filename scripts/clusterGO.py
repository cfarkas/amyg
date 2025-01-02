#!/usr/bin/env python3

"""
Unified script to:
1) Parse annotation + GTF,
2) Build network from GO Jaccard similarity (≥0.2),
3) Detect communities and color top 20,
4) Save gene_network.pdf + clustered_genes.csv,
5) Finally, produce a stacked bar plot per contig with duplication-type counts,
   where bar width is proportional to contig size (deduced from GTF).
   Non-top-20 communities are shown in a light grey background.
   
Note: The original GO-enrichment histogram code and pseudo p-value calculations have been removed.
"""

import sys
import getopt
import os
import re
import random
from collections import defaultdict

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx

from networkx.algorithms.community import greedy_modularity_communities
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import matplotlib.lines as mlines


#############################
# Argument Parsing
#############################

def parse_arguments(argv):
    annotation_file = None
    gtf_file = None
    outdir = "./"

    try:
        opts, args = getopt.getopt(argv, "a:g:o:", ["annotation=", "gtf=", "outdir="])
    except getopt.GetoptError as e:
        print(e)
        print("Usage: clusterGO_fixed.py -a <annotation_file> -g <gtf_file> -o <outdir>")
        sys.exit(2)

    for opt, arg in opts:
        if opt in ("-a", "--annotation"):
            annotation_file = arg
        elif opt in ("-g", "--gtf"):
            gtf_file = arg
        elif opt in ("-o", "--outdir"):
            outdir = arg

    if not annotation_file or not gtf_file:
        print("ERROR: Must specify -a <annotation_file> and -g <gtf_file>.")
        sys.exit(2)

    return {
        "annotation_file": annotation_file,
        "gtf_file": gtf_file,
        "outdir": outdir
    }


#############################
# Minimal GTF Parser
#############################

def parse_gtf_to_dataframe(gtf_file):
    """
    Collect only 'transcript' lines.
    Extract:
       Name = transcript_id
       duplication_type
       start, end, strand, seqname
    Returns a pandas DataFrame.
    """
    rows = []
    with open(gtf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue

            seqname, source, feature, start, end, score, strand, frame, attributes = parts
            if feature != "transcript":
                continue

            attr_dict = {}
            for kv in attributes.split(';'):
                kv = kv.strip()
                if kv:
                    key_val = kv.split(' ', 1)
                    if len(key_val) == 2:
                        key = key_val[0]
                        val = key_val[1].replace('"', '')
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
# Extract GO Terms by Regex
#############################

def extract_go_terms(go_str):
    """
    Find GO terms like GO:1234567 via regex. Return them joined by ';'.
    If none found, return ''.
    """
    matches = re.findall(r'(GO:\d+)', go_str)
    return ";".join(matches)


#############################
# New: Contig Stacked Bar Plot 
#############################

def plot_contig_stacked_bar(merged_df, outdir):
    """
    Creates a stacked bar plot (one bar per contig / 'seqname'):
      - Stacks counts of "ancient", "recent", and "other" duplication
      - Each bar's width is proportional to contig size = max(end) - min(start) + 1
      - Bars are sorted by ascending contig length
    Saves to contig_stacked_bar.pdf in outdir.
    """
    # Group by seqname
    if "seqname" not in merged_df.columns or len(merged_df["seqname"].dropna().unique()) == 0:
        print("No valid seqname data found, skipping contig stacked bar plot.")
        return

    group_data = []
    for seqname, grp in merged_df.groupby("seqname"):
        # contig length
        min_start = grp["start"].min()
        max_end = grp["end"].max()
        contig_len = max_end - min_start + 1

        # duplication counts
        ancient_count = sum(grp["duplication_type"] == "ancient")
        recent_count = sum(grp["duplication_type"] == "recent")
        other_count = len(grp) - (ancient_count + recent_count)

        group_data.append({
            "seqname": seqname,
            "length": contig_len,
            "ancient": ancient_count,
            "recent": recent_count,
            "other": other_count
        })

    if len(group_data) == 0:
        print("No grouped data for contigs, skipping contig stacked bar plot.")
        return

    df_contigs = pd.DataFrame(group_data)
    # sort by length ascending
    df_contigs = df_contigs.sort_values("length")

    # We'll do a stacked bar, with x positions 1..N
    x_positions = np.arange(len(df_contigs))
    bar_widths = df_contigs["length"].values.astype(float)

    # Optional: normalize bar_width so it's not too huge
    max_len = bar_widths.max()
    scale_factor = 1.0 / max_len
    scaled_widths = bar_widths * scale_factor

    # The heights are duplication counts
    ancient_vals = df_contigs["ancient"].values
    recent_vals = df_contigs["recent"].values
    other_vals = df_contigs["other"].values

    plt.figure(figsize=(12, 6))

    # Start with "ancient"
    plt.bar(x_positions, ancient_vals, width=scaled_widths, color="green", label="ancient")

    # Then "recent"
    plt.bar(x_positions, recent_vals, bottom=ancient_vals, width=scaled_widths,
            color="orange", label="recent")

    # Then "other"
    bottom2 = ancient_vals + recent_vals
    plt.bar(x_positions, other_vals, bottom=bottom2, width=scaled_widths,
            color="gray", label="other")

    plt.xticks(x_positions, df_contigs["seqname"], rotation=90)
    plt.xlabel("Contigs (ordered by ascending length)")
    plt.ylabel("Number of Genes")
    plt.legend(loc="upper left")
    plt.title("Stacked Bar: Genes by Duplication Type per Contig (bar width ~ contig size)")

    out_pdf = os.path.join(outdir, "contig_stacked_bar.pdf")
    plt.tight_layout()
    plt.savefig(out_pdf)
    plt.close()
    print(f"Contig stacked bar saved to {out_pdf}")


#############################
# Main Execution
#############################

def main():
    # 1) Parse CLI
    args = parse_arguments(sys.argv[1:])
    annotation_file = args["annotation_file"]
    gtf_file = args["gtf_file"]
    outdir = args["outdir"]

    # 2) Ensure outdir
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    # 3) Read annotation + GTF
    annot_df = pd.read_csv(annotation_file, sep='\t', dtype=str).fillna('')
    gtf_df = parse_gtf_to_dataframe(gtf_file)
    merged_df = pd.merge(annot_df, gtf_df, on="Name", how="left")

    # 4) Extract GO terms
    merged_df['GO_clean'] = merged_df['GO'].fillna('').apply(extract_go_terms)

    # 5) Build Jaccard≥0.2 graph
    G = nx.Graph()
    for i, row in merged_df.iterrows():
        G.add_node(row["Name"], duplication=row["duplication_type"])

    jaccard_threshold = 0.2
    for i in range(len(merged_df)):
        go_i = set(merged_df['GO_clean'].iloc[i].split(';')) if merged_df['GO_clean'].iloc[i] else set()
        for j in range(i+1, len(merged_df)):
            go_j = set(merged_df['GO_clean'].iloc[j].split(';')) if merged_df['GO_clean'].iloc[j] else set()
            union = go_i.union(go_j)
            if not union:
                continue
            inter = go_i.intersection(go_j)
            sim = len(inter)/len(union)
            if sim >= jaccard_threshold:
                G.add_edge(merged_df.iloc[i]["Name"], merged_df.iloc[j]["Name"], weight=sim)

    # 6) Community detection
    if G.number_of_edges() > 0:
        comm_list = list(greedy_modularity_communities(G))
    else:
        comm_list = []

    # map node->community
    node2comm = {}
    for c_id, cset in enumerate(comm_list):
        for n in cset:
            node2comm[n] = c_id

    # for isolated
    for n in G.nodes():
        if n not in node2comm:
            node2comm[n] = -1

    merged_df["community_id"] = merged_df["Name"].map(node2comm).fillna(-1).astype(int)

    # 7) Sort communities by size
    sizes = [(i, len(cset)) for i, cset in enumerate(comm_list)]
    sizes.sort(key=lambda x: x[1], reverse=True)

    # 8) top 20
    top_20 = sizes[:20]
    top20_ids = [t[0] for t in top_20]

    # color map
    cmap20 = plt.cm.get_cmap('tab20', 20)
    comm_color_map = {}
    for rank, (cid, size_) in enumerate(top_20, start=1):
        comm_color_map[cid] = cmap20(rank-1)
    # others => light gray, but make it a bit lighter
    grey_bg = (0.9, 0.9, 0.9, 1.0)  # highlight tab20
    comm_color_map[-1] = grey_bg
    for (cid, csize) in sizes:
        if cid not in top20_ids:
            comm_color_map[cid] = grey_bg

    # 9) Layout
    pos = nx.spring_layout(G, k=0.15, seed=42)

    # 10) subsets for shapes
    tri_x, tri_y, tri_c = [], [], []
    sq_x, sq_y, sq_c = [], [], []
    cir_x, cir_y, cir_c = [], [], []

    for node in G.nodes():
        dup_type = G.nodes[node]["duplication"]
        c_id = node2comm[node]
        color_rgba = comm_color_map.get(c_id, grey_bg)
        x, y = pos[node]
        if dup_type == "ancient":
            tri_x.append(x); tri_y.append(y); tri_c.append(color_rgba)
        elif dup_type == "recent":
            sq_x.append(x); sq_y.append(y); sq_c.append(color_rgba)
        else:
            cir_x.append(x); cir_y.append(y); cir_c.append(color_rgba)

    plt.figure(figsize=(10, 8))
    nx.draw_networkx_edges(G, pos, alpha=0.3)

    plt.scatter(tri_x, tri_y, s=70, marker='^', c=tri_c, edgecolors='k', alpha=0.9)
    plt.scatter(sq_x,  sq_y, s=70, marker='s', c=sq_c, edgecolors='k', alpha=0.9)
    plt.scatter(cir_x, cir_y, s=70, marker='o', c=cir_c, edgecolors='k', alpha=0.9)

    plt.title("Top 20 Communities (colored), shaped by duplication_type")
    plt.axis('off')

    # 11) label top20 with arrow+box
    for rank, (cid, size_) in enumerate(top_20, start=1):
        cset = comm_list[cid] if 0 <= cid < len(comm_list) else []
        if not cset:
            continue
        x_coords = [pos[n][0] for n in cset if n in pos]
        y_coords = [pos[n][1] for n in cset if n in pos]
        if len(x_coords) == 0:
            continue
        xc = sum(x_coords)/len(x_coords)
        yc = sum(y_coords)/len(y_coords)
        x_text, y_text = xc+0.03, yc+0.03
        plt.annotate(str(rank),
                     xy=(xc,yc),
                     xytext=(x_text,y_text),
                     arrowprops=dict(arrowstyle="->", color="black", lw=1),
                     bbox=dict(boxstyle="square,pad=0.25", fc="white", ec="black", alpha=0.9),
                     fontsize=11, ha='center', va='center', color='black')

    # 12) legend
    shape_anc = mlines.Line2D([], [], color='black', marker='^', linestyle='None',
                              markersize=8, label='ancient duplication')
    shape_rec = mlines.Line2D([], [], color='black', marker='s', linestyle='None',
                              markersize=8, label='recent duplication')
    shape_oth = mlines.Line2D([], [], color='black', marker='o', linestyle='None',
                              markersize=8, label='no_duplication/other')
    shape_handles = [shape_anc, shape_rec, shape_oth]

    comm_handles = []
    for rank, (cid, size_) in enumerate(top_20, start=1):
        rgba = comm_color_map.get(cid, grey_bg)
        lbl = f"{rank} (n={size_})"
        patch = mpatches.Patch(facecolor=rgba, edgecolor='none', label=lbl)
        comm_handles.append(patch)

    all_handles = shape_handles + comm_handles

    plt.legend(handles=all_handles,
               loc='upper center',
               bbox_to_anchor=(0.5, -0.15),
               fancybox=True,
               shadow=True,
               ncol=3)

    # 13) Save network
    pdf_path = os.path.join(outdir, "gene_network.pdf")
    plt.savefig(pdf_path, bbox_inches='tight')
    plt.close()

    csv_path = os.path.join(outdir, "clustered_genes.csv")
    merged_df.to_csv(csv_path, index=False)
    print(f"Done! PDF -> {pdf_path}\nCSV -> {csv_path}\n")

    # 14) produce contig stacked bar figure
    plot_contig_stacked_bar(merged_df, outdir)


if __name__ == "__main__":
    main()
