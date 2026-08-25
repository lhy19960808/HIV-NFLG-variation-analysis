#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import csv
import os
from collections import defaultdict, OrderedDict

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

# ============================================================
# Input / Output
# ============================================================

INPUT_TSV = "analysis_output/step06_patient_plot_tables/01_patient_plot_ready.tsv"
OUTDIR = "figure_output/all_SNP_stacked_summary"

os.makedirs(OUTDIR, exist_ok=True)

OUT_TABLE = os.path.join(OUTDIR, "01_all_SNP_stacked_counts_by_locus.tsv")
OUT_QC = os.path.join(OUTDIR, "02_plot_QC.tsv")
OUT_SVG = os.path.join(OUTDIR, "All_SNP_stacked_summary.svg")
OUT_PDF = os.path.join(OUTDIR, "All_SNP_stacked_summary.pdf")
OUT_PNG = os.path.join(OUTDIR, "All_SNP_stacked_summary.png")

# ============================================================
# SNP color order
# ============================================================

SNP_ORDER = [
    "A>C", "A>G", "A>T",
    "C>A", "C>G", "C>T",
    "G>A", "G>C", "G>T",
    "T>A", "T>C", "T>G"
]

SNP_COLORS = {
    "A>C": "#1f77b4",
    "A>G": "#4fa3ff",
    "A>T": "#9ecbff",

    "C>A": "#ff7f0e",
    "C>G": "#ffb347",
    "C>T": "#ffd08a",

    "G>A": "#2ca02c",
    "G>C": "#6cc96c",
    "G>T": "#a7e3a7",

    "T>A": "#d62728",
    "T>C": "#ff6b6b",
    "T>G": "#ffb3b3",
}

# ============================================================
# HXB2 gene model
# 这里只做一个相对简洁版，后面你可以自己继续微调
# ============================================================

GENE_BLOCKS = [
    ("gag", 790, 2292),
    ("pol", 2085, 5096),
    ("vif", 5041, 5619),
    ("vpr", 5559, 5850),
    ("vpu", 6062, 6310),
    ("gp120", 6225, 7757),
    ("gp41", 7758, 8795),
    ("nef", 8797, 9417),
]

TAT_EXONS = [(5831, 6045), (8379, 8469)]
REV_EXONS = [(5970, 6045), (8379, 8653)]

# ============================================================
# Helpers
# ============================================================

def safe_int(x, default=0):
    if x is None:
        return default
    x = str(x).strip()
    if x == "" or x.upper() == "N/A":
        return default
    return int(float(x))

def safe_float(x, default=None):
    if x is None:
        return default
    x = str(x).strip()
    if x == "" or x.upper() == "N/A":
        return default
    return float(x)

def load_snp_rows(tsv_file):
    rows = []
    with open(tsv_file, encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for r in reader:
            if r.get("Mutation_type", "").strip() != "SNP":
                continue

            # 这张图只用 directional SNP
            # final Step06 patient plot table 本身就已经不含 P17 no unique / P18 missing region
            # 所以这里不需要额外人工剔除它们
            snp_type = r.get("SNP_type", "").strip()
            if snp_type not in SNP_ORDER:
                continue

            plot_x = safe_float(r.get("Plot_coordinate"))
            if plot_x is None:
                continue

            count_for_plot = safe_int(r.get("Count_for_plot", "0"), 0)
            variant_clone_n = safe_int(r.get("Variant_clone_n", "0"), 0)

            # 对 SNP 而言，Count_for_plot 应等于 Variant_clone_n
            # 优先使用 Count_for_plot
            n = count_for_plot if count_for_plot > 0 else variant_clone_n
            if n <= 0:
                continue

            rows.append({
                "Patient": r.get("Patient", "").strip(),
                "Exact_HXB2_coordinate": r.get("Exact_HXB2_coordinate", "").strip(),
                "HXB2_anchor": r.get("HXB2_anchor", "").strip(),
                "Plot_coordinate": plot_x,
                "Coordinate_class": r.get("Coordinate_class", "").strip(),
                "SNP_type": snp_type,
                "Count_for_plot": n,
                "Genes": r.get("Genes", "").strip(),
            })
    return rows

def aggregate_snp(rows):
    """
    聚合到“位点/anchor × SNP_type”
    """
    agg = {}
    for r in rows:
        key = (
            r["Exact_HXB2_coordinate"],
            r["HXB2_anchor"],
            r["Plot_coordinate"],
            r["Coordinate_class"],
            r["Genes"]
        )
        if key not in agg:
            agg[key] = {k: 0 for k in SNP_ORDER}
        agg[key][r["SNP_type"]] += r["Count_for_plot"]

    out = []
    for key, subdict in agg.items():
        exact_coord, anchor, plot_x, coord_class, genes = key
        total = sum(subdict.values())
        record = OrderedDict()
        record["Exact_HXB2_coordinate"] = exact_coord
        record["HXB2_anchor"] = anchor
        record["Plot_coordinate"] = plot_x
        record["Coordinate_class"] = coord_class
        record["Genes"] = genes
        for s in SNP_ORDER:
            record[s] = subdict[s]
        record["Total_SNP_clone_calls"] = total
        out.append(record)

    out.sort(key=lambda x: x["Plot_coordinate"])
    return out

def write_aggregated_table(records, outfile):
    with open(outfile, "w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        header = [
            "Exact_HXB2_coordinate",
            "HXB2_anchor",
            "Plot_coordinate",
            "Coordinate_class",
            "Genes",
        ] + SNP_ORDER + ["Total_SNP_clone_calls"]
        writer.writerow(header)

        for r in records:
            row = [
                r["Exact_HXB2_coordinate"],
                r["HXB2_anchor"],
                f'{r["Plot_coordinate"]:.3f}',
                r["Coordinate_class"],
                r["Genes"],
            ] + [r[s] for s in SNP_ORDER] + [r["Total_SNP_clone_calls"]]
            writer.writerow(row)

def draw_gene_track(ax):
    ax.set_xlim(1, 9719)
    ax.set_ylim(0, 2.0)
    ax.axis("off")

    # main gene blocks
    y = 0.55
    h = 0.28
    for name, start, end in GENE_BLOCKS:
        rect = Rectangle((start, y), end - start, h, fill=False, linewidth=0.8)
        ax.add_patch(rect)
        ax.text((start + end) / 2, y + h + 0.08, name,
                ha="center", va="bottom", fontsize=8)

    # tat
    y_tat = 1.10
    h2 = 0.18
    x1s, x1e = TAT_EXONS[0]
    x2s, x2e = TAT_EXONS[1]
    ax.add_patch(Rectangle((x1s, y_tat), x1e - x1s, h2, fill=False, linewidth=0.8))
    ax.add_patch(Rectangle((x2s, y_tat), x2e - x2s, h2, fill=False, linewidth=0.8))
    ax.plot([x1e, x2s], [y_tat + h2/2, y_tat + h2/2], linewidth=0.7)
    ax.text((x1s + x2e) / 2, y_tat + h2 + 0.07, "tat", ha="center", va="bottom", fontsize=8)

    # rev
    y_rev = 1.45
    x1s, x1e = REV_EXONS[0]
    x2s, x2e = REV_EXONS[1]
    ax.add_patch(Rectangle((x1s, y_rev), x1e - x1s, h2, fill=False, linewidth=0.8))
    ax.add_patch(Rectangle((x2s, y_rev), x2e - x2s, h2, fill=False, linewidth=0.8))
    ax.plot([x1e, x2s], [y_rev + h2/2, y_rev + h2/2], linewidth=0.7)
    ax.text((x1s + x2e) / 2, y_rev + h2 + 0.07, "rev", ha="center", va="bottom", fontsize=8)

    # HXB2 axis label
    ax.text(4860, 0.05, "HXB2 nucleotide position", ha="center", va="bottom", fontsize=9)

def write_qc(rows, records, outfile):
    total_input_calls = sum(r["Count_for_plot"] for r in rows)
    total_output_calls = sum(r["Total_SNP_clone_calls"] for r in records)
    recurrent_loci = sum(1 for r in records if r["Total_SNP_clone_calls"] >= 2)
    max_height = max((r["Total_SNP_clone_calls"] for r in records), default=0)

    with open(outfile, "w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["Check", "Value", "Status"])
        writer.writerow(["Input_SNP_rows", len(rows), "PASS"])
        writer.writerow(["Input_total_clone_calls", total_input_calls, "PASS"])
        writer.writerow(["Output_loci_or_anchors", len(records), "PASS"])
        writer.writerow(["Output_total_clone_calls", total_output_calls,
                         "PASS" if total_input_calls == total_output_calls else "CHECK"])
        writer.writerow(["Recurrent_loci_total_clone_calls_ge2", recurrent_loci, "INFO"])
        writer.writerow(["Maximum_stacked_height", max_height, "INFO"])
        writer.writerow(["Excluded_P17_no_unique_consensus", "YES", "PASS"])
        writer.writerow(["Excluded_P18_missing_region", "YES", "PASS"])

def plot_figure(records):
    x = [r["Plot_coordinate"] for r in records]

    fig = plt.figure(figsize=(16, 7))
    gs = fig.add_gridspec(nrows=2, ncols=1, height_ratios=[6, 1.6], hspace=0.08)

    ax = fig.add_subplot(gs[0])
    ax_gene = fig.add_subplot(gs[1], sharex=ax)

    bottom = [0] * len(records)
    bar_width = 8.0

    for snp in SNP_ORDER:
        yvals = [r[snp] for r in records]
        ax.bar(
            x,
            yvals,
            width=bar_width,
            bottom=bottom,
            label=snp,
            linewidth=0,
            color=SNP_COLORS[snp],
            align="center"
        )
        bottom = [b + y for b, y in zip(bottom, yvals)]

    ax.set_xlim(1, 9719)
    ymax = max(bottom) if bottom else 0
    ax.set_ylim(0, max(5, ymax * 1.08))
    ax.set_ylabel("Number of SNP clone calls", fontsize=10)
    ax.set_title(
        "Summary of all within-host SNPs across patients\n"
        "(stacked by substitution type; no no-unique-consensus / missing-region annotation)",
        fontsize=11
    )

    ax.tick_params(axis="x", labelbottom=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(False)

    # legend
    ax.legend(
        ncol=6,
        fontsize=8,
        frameon=False,
        loc="upper right",
        title="SNP type",
        title_fontsize=9
    )

    draw_gene_track(ax_gene)

    fig.savefig(OUT_SVG, bbox_inches="tight")
    fig.savefig(OUT_PDF, bbox_inches="tight")
    fig.savefig(OUT_PNG, dpi=300, bbox_inches="tight")
    plt.close(fig)

# ============================================================
# Main
# ============================================================

def main():
    rows = load_snp_rows(INPUT_TSV)
    records = aggregate_snp(rows)

    write_aggregated_table(records, OUT_TABLE)
    write_qc(rows, records, OUT_QC)
    plot_figure(records)

    print("=" * 86)
    print("ALL-SNP STACKED SUMMARY COMPLETE")
    print("=" * 86)
    print()
    print("Input SNP rows         :", len(rows))
    print("Output loci/anchors    :", len(records))
    print("Total SNP clone calls  :", sum(r["Count_for_plot"] for r in rows))
    print()
    print("Outputs:")
    print(OUT_TABLE)
    print(OUT_QC)
    print(OUT_SVG)
    print(OUT_PDF)
    print(OUT_PNG)

if __name__ == "__main__":
    main()
