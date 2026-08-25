#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Figure 2
Cross-patient distribution of within-host HIV-1 variable sites

Plotting design:
- No biological/statistical data changes.
- SNP stems use one consistent style.
- Recurrent SNP sites (Patient_n >= 2) use one consistent point style.
- DEL markers are unified within one marker track.
- INS markers are unified within one marker track.
- tat/rev connector lines are unified.
- Shorter title.
- Removed the long explanatory sentence from inside the figure.
- Compressed HXB2 gene model vertically.
- HXB2 6403 annotation simplified to "HXB2 6403 / n = 3".

Figure 2 includes ONLY:
  SNP
  DEL
  INS

Excluded:
  No unique consensus
  P17 unresolved structural region
  P18 large missing region
  terminal no-callable regions

Primary cross-patient metric:
  Patient_n
"""

from pathlib import Path
import csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D

ROOT = Path(".")
DATA = (
    ROOT
    / "analysis_output/step08_figure2_plot_data"
    / "01_figure2_plotting_summary.tsv"
)
QC = (
    ROOT
    / "analysis_output/step08_figure2_plot_data"
    / "05_figure2_plot_QC.tsv"
)

OUTDIR = ROOT / "figure_output/figure2"
OUTDIR.mkdir(parents=True, exist_ok=True)

OUT_SVG = OUTDIR / "Figure2_cross_patient_variable_sites.svg"
OUT_PDF = OUTDIR / "Figure2_cross_patient_variable_sites.pdf"
OUT_PNG = OUTDIR / "Figure2_cross_patient_variable_sites.png"
OUT_QC = OUTDIR / "Figure2_plot_QC.tsv"

HXB2_MIN = 1
HXB2_MAX = 9719

# Exact HXB2 gene coordinates used throughout the project.
GENES = [
    ("gag",   790, 2292, 0),
    ("pol",  2085, 5096, 1),
    ("vif",  5041, 5619, 0),
    ("vpr",  5559, 5850, 1),
    ("vpu",  6062, 6310, 0),
    ("gp120",6225, 7757, 1),
    ("gp41", 7758, 8795, 0),
    ("nef",  8797, 9417, 1),
]

SPLIT_GENES = {
    "tat": [(5831, 6045), (8379, 8469)],
    "rev": [(5970, 6045), (8379, 8653)],
}


def read_tsv(path):
    if not path.exists():
        raise FileNotFoundError(f"Missing required file: {path}")
    with open(path, encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def write_tsv(path, fields, rows):
    with open(path, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(
            f,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n"
        )
        w.writeheader()
        for r in rows:
            w.writerow(r)


# ============================================================
# Validate final Step08 Figure 2 plotting data
# ============================================================

qc_rows = read_tsv(QC)

ready = [
    r for r in qc_rows
    if r["Check"] == "READY_FOR_FIGURE2_PLOTTING"
]

if not ready or ready[0]["Status"] != "PASS":
    raise RuntimeError(
        "Step08 plot data are not READY_FOR_FIGURE2_PLOTTING. Check Step08 QC first."
    )

rows = read_tsv(DATA)

bad = [
    r for r in rows
    if r["Mutation_type"] not in {"SNP", "DEL", "INS"}
]

if bad:
    raise RuntimeError(
        "Non-SNP/DEL/INS records detected in Figure 2 input."
    )

snps = [r for r in rows if r["Mutation_type"] == "SNP"]
dels = [r for r in rows if r["Mutation_type"] == "DEL"]
ins = [r for r in rows if r["Mutation_type"] == "INS"]

if len(snps) != 973:
    raise RuntimeError(
        f"Expected 973 SNP loci/anchors, observed {len(snps)}"
    )

if len(dels) != 32:
    raise RuntimeError(
        f"Expected 32 DEL loci, observed {len(dels)}"
    )

if len(ins) != 13:
    raise RuntimeError(
        f"Expected 13 INS anchors, observed {len(ins)}"
    )


# ============================================================
# Prepare plot arrays
# ============================================================

snp_x = [float(r["Plot_coordinate"]) for r in snps]
snp_y = [int(r["Patient_n"]) for r in snps]

recurrent = [
    r for r in snps
    if int(r["Patient_n"]) >= 2
]

recurrent_x = [
    float(r["Plot_coordinate"])
    for r in recurrent
]

recurrent_y = [
    int(r["Patient_n"])
    for r in recurrent
]

del_x = [
    float(r["Plot_coordinate"])
    for r in dels
]

ins_x = [
    float(r["Plot_coordinate"])
    for r in ins
]


# ============================================================
# Figure layout
# ============================================================

fig = plt.figure(figsize=(14.2, 5.7))

# Main mutation panel
ax = fig.add_axes([0.075, 0.31, 0.895, 0.51])

# Compact HXB2 gene model
gene_ax = fig.add_axes([0.075, 0.105, 0.895, 0.145], sharex=ax)


# ============================================================
# SNP track
# ============================================================

# One vlines call keeps all SNP stems visually consistent.
# Height remains exactly Patient_n.
line_widths = [
    0.48 if y == 1 else
    1.05 if y == 2 else
    1.55
    for y in snp_y
]

ax.vlines(
    snp_x,
    0,
    snp_y,
    linewidths=line_widths,
    alpha=0.72,
    zorder=2
)

# One scatter call gives every recurrent site the same marker style.
ax.scatter(
    recurrent_x,
    recurrent_y,
    marker="o",
    s=12,
    zorder=4
)


# ============================================================
# DEL / INS marker tracks
# ============================================================

DEL_Y = 4.00
INS_Y = 4.62

# Single calls per type prevent automatic within-type color cycling.
ax.scatter(
    del_x,
    [DEL_Y] * len(del_x),
    marker="D",
    s=18,
    zorder=5
)

ax.scatter(
    ins_x,
    [INS_Y] * len(ins_x),
    marker="v",
    s=24,
    zorder=5
)


# ============================================================
# HXB2 6403 annotation
# ============================================================

three_patient_sites = [
    r for r in snps
    if int(r["Patient_n"]) == 3
]

for r in three_patient_sites:
    x = float(r["Plot_coordinate"])

    ax.annotate(
        f'HXB2 {r["HXB2_locus_or_anchor"]}\nn = 3',
        xy=(x, 3),
        xytext=(x + 170, 3.32),
        fontsize=7.5,
        ha="left",
        va="bottom",
        arrowprops={
            "arrowstyle": "-",
            "linewidth": 0.55
        }
    )


# ============================================================
# Main axes formatting
# ============================================================

ax.set_xlim(HXB2_MIN, HXB2_MAX)
ax.set_ylim(0, 5.08)

ax.set_ylabel(
    "Patient count",
    fontsize=10.5
)

ax.set_yticks(
    [0, 1, 2, 3, DEL_Y, INS_Y]
)

ax.set_yticklabels(
    ["0", "1", "2", "3", "DEL", "INS"]
)

xticks = [
    1, 1000, 2000, 3000, 4000,
    5000, 6000, 7000, 8000, 9000, 9719
]

ax.set_xticks(xticks)
ax.tick_params(
    axis="x",
    labelbottom=False
)

ax.tick_params(
    axis="both",
    width=0.55,
    length=3,
    labelsize=8.5
)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["left"].set_linewidth(0.55)
ax.spines["bottom"].set_linewidth(0.55)


# ============================================================
# Figure legend
# ============================================================

legend_handles = [
    Line2D(
        [0], [0],
        marker="|",
        linestyle="None",
        markersize=12,
        label="SNP locus"
    ),
    Line2D(
        [0], [0],
        marker="o",
        linestyle="None",
        markersize=4,
        label="Recurrent SNP locus (Patient_n ≥ 2)"
    ),
    Line2D(
        [0], [0],
        marker="D",
        linestyle="None",
        markersize=5,
        label="DEL locus"
    ),
    Line2D(
        [0], [0],
        marker="v",
        linestyle="None",
        markersize=5,
        label="INS anchor"
    ),
]

ax.legend(
    handles=legend_handles,
    ncol=4,
    loc="upper center",
    bbox_to_anchor=(0.5, 1.14),
    frameon=False,
    fontsize=8,
    handletextpad=0.45,
    columnspacing=1.15
)


# ============================================================
# HXB2 gene model
# ============================================================

gene_ax.set_xlim(HXB2_MIN, HXB2_MAX)
gene_ax.set_ylim(0, 2.42)

gene_h = 0.31

primary_track_y = {
    0: 0.28,
    1: 0.76
}

for name, start, end, track in GENES:
    y = primary_track_y[track]

    gene_ax.add_patch(
        Rectangle(
            (start, y),
            end - start + 1,
            gene_h,
            fill=False,
            linewidth=0.58
        )
    )

    gene_ax.text(
        (start + end) / 2,
        y + gene_h / 2,
        name,
        ha="center",
        va="center",
        fontsize=7.2
    )


# Split tat/rev exons
split_y = {
    "tat": 1.28,
    "rev": 1.78
}

connector_y = []
connector_start = []
connector_end = []

for gene, exons in SPLIT_GENES.items():
    y = split_y[gene]
    (s1, e1), (s2, e2) = exons

    gene_ax.add_patch(
        Rectangle(
            (s1, y),
            e1 - s1 + 1,
            0.25,
            fill=False,
            linewidth=0.58
        )
    )

    gene_ax.add_patch(
        Rectangle(
            (s2, y),
            e2 - s2 + 1,
            0.25,
            fill=False,
            linewidth=0.58
        )
    )

    connector_y.append(y + 0.125)
    connector_start.append(e1)
    connector_end.append(s2)

    gene_ax.text(
        (e1 + s2) / 2,
        y + 0.30,
        gene,
        ha="center",
        va="bottom",
        fontsize=7.2
    )

# One hlines call keeps tat/rev connectors visually consistent.
gene_ax.hlines(
    connector_y,
    connector_start,
    connector_end,
    linewidth=0.50
)


# ============================================================
# HXB2 axis
# ============================================================

gene_ax.set_xlabel(
    "HXB2 nucleotide position",
    fontsize=10
)

gene_ax.set_xticks(xticks)

gene_ax.tick_params(
    axis="x",
    width=0.55,
    length=3,
    labelsize=8.5
)

gene_ax.set_yticks([])

gene_ax.spines["top"].set_visible(False)
gene_ax.spines["right"].set_visible(False)
gene_ax.spines["left"].set_visible(False)
gene_ax.spines["bottom"].set_linewidth(0.55)


# ============================================================
# Title
# ============================================================

fig.suptitle(
    "Cross-patient distribution of within-host HIV-1 variable sites",
    y=0.965,
    fontsize=13
)


# ============================================================
# Save
# ============================================================

fig.savefig(
    OUT_SVG,
    bbox_inches="tight"
)

fig.savefig(
    OUT_PDF,
    bbox_inches="tight"
)

fig.savefig(
    OUT_PNG,
    dpi=600,
    bbox_inches="tight"
)

plt.close(fig)


# ============================================================
# QC
# ============================================================

max_snp_patient_n = max(
    int(r["Patient_n"])
    for r in snps
)

recurrent_n = sum(
    int(r["Patient_n"]) >= 2
    for r in snps
)

same_substitution_n = sum(
    r["Substitution_pattern"] == "SAME_SUBSTITUTION"
    for r in recurrent
)

mixed_substitution_n = sum(
    r["Substitution_pattern"] == "MIXED_SUBSTITUTIONS"
    for r in recurrent
)

qc_out = [
    {
        "Check": "SNP_locus_or_anchor_n",
        "Observed": len(snps),
        "Expected": 973,
        "Status": "PASS" if len(snps) == 973 else "FAIL",
    },
    {
        "Check": "Recurrent_SNP_locus_n",
        "Observed": recurrent_n,
        "Expected": 52,
        "Status": "PASS" if recurrent_n == 52 else "FAIL",
    },
    {
        "Check": "Same_substitution_recurrent_n",
        "Observed": same_substitution_n,
        "Expected": 21,
        "Status": "PASS" if same_substitution_n == 21 else "FAIL",
    },
    {
        "Check": "Mixed_substitution_recurrent_n",
        "Observed": mixed_substitution_n,
        "Expected": 31,
        "Status": "PASS" if mixed_substitution_n == 31 else "FAIL",
    },
    {
        "Check": "Maximum_SNP_Patient_n",
        "Observed": max_snp_patient_n,
        "Expected": 3,
        "Status": "PASS" if max_snp_patient_n == 3 else "FAIL",
    },
    {
        "Check": "DEL_locus_n",
        "Observed": len(dels),
        "Expected": 32,
        "Status": "PASS" if len(dels) == 32 else "FAIL",
    },
    {
        "Check": "INS_anchor_n",
        "Observed": len(ins),
        "Expected": 13,
        "Status": "PASS" if len(ins) == 13 else "FAIL",
    },
    {
        "Check": "Only_SNP_DEL_INS",
        "Observed": len(bad),
        "Expected": 0,
        "Status": "PASS" if len(bad) == 0 else "FAIL",
    },
    {
        "Check": "Three_patient_SNP_locus_n",
        "Observed": len(three_patient_sites),
        "Expected": 1,
        "Status": "PASS" if len(three_patient_sites) == 1 else "FAIL",
    },
    {
        "Check": "SVG_created",
        "Observed": str(OUT_SVG),
        "Expected": "created",
        "Status": "PASS" if OUT_SVG.exists() else "FAIL",
    },
    {
        "Check": "PDF_created",
        "Observed": str(OUT_PDF),
        "Expected": "created",
        "Status": "PASS" if OUT_PDF.exists() else "FAIL",
    },
    {
        "Check": "PNG_created",
        "Observed": str(OUT_PNG),
        "Expected": "created",
        "Status": "PASS" if OUT_PNG.exists() else "FAIL",
    },
]

write_tsv(
    OUT_QC,
    ["Check", "Observed", "Expected", "Status"],
    qc_out
)


# ============================================================
# Console summary
# ============================================================

print("=" * 84)
print("FIGURE 2 COMPLETE")
print("=" * 84)

print("SNP loci/anchors              :", len(snps))
print("Recurrent SNP loci >=2        :", recurrent_n)
print("Same-substitution recurrent   :", same_substitution_n)
print("Mixed-substitution recurrent  :", mixed_substitution_n)
print("Maximum SNP Patient_n         :", max_snp_patient_n)
print("DEL loci                      :", len(dels))
print("INS anchors                   :", len(ins))

print()
print("SVG:", OUT_SVG)
print("PDF:", OUT_PDF)
print("PNG:", OUT_PNG)
print("QC :", OUT_QC)
