#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import csv
import re
from collections import defaultdict, OrderedDict
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Polygon, Rectangle

DATA_DIR = Path("analysis_output") / "step06_patient_plot_tables"
PLOT_TABLE = DATA_DIR / "01_patient_plot_ready.tsv"
META_TABLE = DATA_DIR / "05_patient_panel_metadata.tsv"
NOUNIQUE_TABLE = DATA_DIR / "06_no_unique_consensus_regions.tsv"
PLOT_QC = DATA_DIR / "04_plot_data_QC.tsv"

OUTDIR = Path("figure_output") / "figure1"
OUTDIR.mkdir(parents=True, exist_ok=True)

OUT_SVG = OUTDIR / "Figure1_patient_mutation_landscape.svg"
OUT_PDF = OUTDIR / "Figure1_patient_mutation_landscape.pdf"
OUT_PNG = OUTDIR / "Figure1_patient_mutation_landscape.png"
OUT_QC = OUTDIR / "Figure1_plot_QC.tsv"

XMIN, XMAX = 1, 9719
YMIN, YMAX = -0.8, 25.0

DEL_Y = 18.2
INS_Y = 21.0
DEL_LABEL_LEVELS = [19.15, 20.05]
INS_LABEL_LEVELS = [21.95, 22.85, 23.75]
EVENT_LABEL_MIN_X_SEPARATION = 150.0

NOUNIQUE_Y = 8.3
UNRESOLVED_Y = 5.8
UNRESOLVED_HEIGHT = 0.82

BAR_WIDTH = 0.55
DEL_MARKER_SIZE = 13
INS_MARKER_SIZE = 15
NOUNIQUE_MARKER_SIZE = 14

EVENT_LABEL_FS = 5.0
PATIENT_LABEL_FS = 8.8
RANGE_LABEL_FS = 6.2
TICK_LABEL_FS = 6.2
GENE_LABEL_FS = 5.8
TITLE_FS = 11.5

AXIS_COLOR = "#9A9A9A"
AXIS_TEXT_COLOR = "#666666"
BASELINE_LW = 0.34
TICK_LW = 0.30

DEL_EDGE = "#222222"
INS_EDGE = "#222222"
NOUNIQUE_COLOR = "#111111"
UNRESOLVED_COLOR = "#111111"

SNP_ORDER = [
    "A>C", "A>G", "A>T",
    "C>A", "C>G", "C>T",
    "G>A", "G>C", "G>T",
    "T>A", "T>C", "T>G",
]

SNP_COLORS = OrderedDict({
    "A>C": "#1f77b4",
    "A>G": "#2ca02c",
    "A>T": "#17becf",
    "C>A": "#ff7f0e",
    "C>G": "#bcbd22",
    "C>T": "#d62728",
    "G>A": "#9467bd",
    "G>C": "#8c564b",
    "G>T": "#e377c2",
    "T>A": "#7f7f7f",
    "T>C": "#aec7e8",
    "T>G": "#98df8a",
})

PRIMARY_TRACK_1 = [
    ("gag", 790, 2292),
    ("vif", 5041, 5619),
    ("vpu", 6062, 6310),
    ("gp41", 7758, 8795),
]

PRIMARY_TRACK_2 = [
    ("pol", 2085, 5096),
    ("vpr", 5559, 5850),
    ("gp120", 6225, 7757),
    ("nef", 8797, 9417),
]

SPLIT_GENES = {
    "tat": [(5831, 6045), (8379, 8469)],
    "rev": [(5970, 6045), (8379, 8653)],
}

def require(path):
    if not path.exists():
        raise FileNotFoundError(path)

def read_tsv(path):
    with open(path, "r", encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))

def write_tsv(path, fields, rows):
    with open(path, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for row in rows:
            w.writerow({k: row.get(k, "") for k in fields})

def safe_int(x, default=None):
    try:
        s = str(x).strip()
        return default if s == "" else int(float(s))
    except Exception:
        return default

def safe_float(x, default=None):
    try:
        s = str(x).strip()
        return default if s == "" else float(s)
    except Exception:
        return default

def extract_numbers(text):
    return [int(x) for x in re.findall(r"\d+", str(text or ""))]

def interval_from_anchors(left_text, right_text):
    left_nums = extract_numbers(left_text)
    right_nums = extract_numbers(right_text)
    left = left_nums[0] if left_nums else None
    right = right_nums[-1] if right_nums else None
    if left is None and right is None:
        return None, None
    if left is None:
        left = right
    if right is None:
        right = left
    if left > right:
        left, right = right, left
    return float(left), float(right)

def patient_sort_key(patient):
    m = re.fullmatch(r"CN(\d{4})AH(\d+)-(\d+)", str(patient))
    return tuple(map(int, m.groups())) if m else (9999, 9999, str(patient))

def assign_label_levels(event_rows, levels, min_sep):
    items = []
    for r in event_rows:
        x = safe_float(r.get("Plot_coordinate"))
        count = safe_int(r.get("Count_for_plot"), 0)
        if x is not None and count and count > 0:
            items.append({"row": r, "x": x, "count": count})
    items.sort(key=lambda z: z["x"])
    last_x = [-10**12 for _ in levels]
    for item in items:
        x = item["x"]
        chosen = None
        for i in range(len(levels)):
            if x - last_x[i] >= min_sep:
                chosen = i
                break
        if chosen is None:
            chosen = max(range(len(levels)), key=lambda i: x - last_x[i])
        item["label_y"] = levels[chosen]
        last_x[chosen] = x
    return items

def style_patient_axis(ax):
    ax.set_xlim(XMIN, XMAX)
    ax.set_ylim(YMIN, YMAX)
    for side in ["top", "right", "left", "bottom"]:
        ax.spines[side].set_visible(False)
    ax.axhline(0, color=AXIS_COLOR, lw=BASELINE_LW, zorder=1)
    ax.set_yticks([0, 5, 10, 15, 20])
    ax.set_yticklabels([0, 5, 10, 15, 20], fontsize=TICK_LABEL_FS, color=AXIS_TEXT_COLOR)
    ax.tick_params(axis="y", width=TICK_LW, length=1.8, color=AXIS_COLOR, pad=1.8)
    ax.set_xticks([])
    ax.tick_params(axis="x", length=0)

def draw_patient_label(ax, patient, start, end):
    ax.text(-0.055, 0.73, patient, transform=ax.transAxes,
            ha="right", va="center", fontsize=PATIENT_LABEL_FS,
            fontweight="bold", clip_on=False)
    ax.text(-0.055, 0.545, f"HXB2 {start}-{end}", transform=ax.transAxes,
            ha="right", va="center", fontsize=RANGE_LABEL_FS,
            color=AXIS_TEXT_COLOR, clip_on=False)

def draw_snp_bars(ax, rows):
    grouped = defaultdict(lambda: defaultdict(int))
    for r in rows:
        x = safe_float(r.get("Plot_coordinate"))
        snp_type = str(r.get("SNP_type", "")).strip()
        count = safe_int(r.get("Count_for_plot"), 0)
        if x is not None and snp_type and count > 0:
            grouped[x][snp_type] += count
    max_height = 0
    for x in sorted(grouped):
        bottom = 0
        for snp_type in SNP_ORDER:
            h = grouped[x].get(snp_type, 0)
            if h <= 0:
                continue
            ax.bar(x, h, width=BAR_WIDTH, bottom=bottom,
                   color=SNP_COLORS[snp_type], edgecolor="none",
                   linewidth=0, align="center", zorder=3)
            bottom += h
        max_height = max(max_height, bottom)
    return max_height

def draw_event_track(ax, rows, marker_y, marker, size, label_levels):
    placed = assign_label_levels(rows, label_levels, EVENT_LABEL_MIN_X_SEPARATION)
    for item in placed:
        x, count, label_y = item["x"], item["count"], item["label_y"]
        edge = DEL_EDGE if marker == "D" else INS_EDGE
        ax.scatter([x], [marker_y], marker=marker, s=size,
                   facecolors="white", edgecolors=edge,
                   linewidths=0.55, zorder=5)
        ax.plot([x, x], [marker_y + 0.18, label_y - 0.16],
                color="#777777", lw=0.30, zorder=4)
        ax.text(x, label_y, str(count), ha="center", va="bottom",
                fontsize=EVENT_LABEL_FS, color="#222222",
                bbox=dict(facecolor="white", edgecolor="none", alpha=0.90, pad=0.10),
                zorder=6)
    return placed

def draw_hex_band(ax, start, end, y, height=UNRESOLVED_HEIGHT):
    if start is None or end is None:
        return
    if end < start:
        start, end = end, start
    width = end - start
    if width <= 8:
        ax.scatter([(start + end) / 2.0], [y], marker="D",
                   s=NOUNIQUE_MARKER_SIZE, color=NOUNIQUE_COLOR, zorder=5)
        return
    taper = min(60.0, max(10.0, width * 0.035))
    taper = min(taper, width * 0.30)
    pts = [
        (start, y),
        (start + taper, y + height / 2.0),
        (end - taper, y + height / 2.0),
        (end, y),
        (end - taper, y - height / 2.0),
        (start + taper, y - height / 2.0),
    ]
    ax.add_patch(Polygon(
        pts, closed=True, facecolor=UNRESOLVED_COLOR,
        edgecolor=UNRESOLVED_COLOR, linewidth=0.35,
        alpha=0.92, zorder=4
    ))

def draw_no_unique(ax, rows):
    point_n = 0
    band_n = 0
    for r in rows:
        cls = str(r.get("No_unique_class", "")).strip()
        length = safe_int(r.get("Length_alignment_columns"), 1)
        start, end = interval_from_anchors(
            r.get("Left_anchor_HXB2_coordinate", ""),
            r.get("Right_anchor_HXB2_coordinate", "")
        )
        if start is None or end is None:
            continue
        if "base_only" in cls.lower() and length is not None and length <= 3:
            x = (start + end) / 2.0
            ax.scatter([x], [NOUNIQUE_Y], marker="D",
                       s=NOUNIQUE_MARKER_SIZE, color=NOUNIQUE_COLOR,
                       linewidths=0, zorder=5)
            point_n += 1
        else:
            draw_hex_band(ax, start, end, UNRESOLVED_Y)
            band_n += 1

    if point_n > 0:
        ax.text(0.01, 0.355, "No unique consensus",
                transform=ax.transAxes, ha="left", va="center",
                fontsize=5.8, color="#222222")
    if band_n > 0:
        ax.text(0.01, 0.255, "Unresolved region",
                transform=ax.transAxes, ha="left", va="center",
                fontsize=5.6, color="#222222")
    return {"point_n": point_n, "band_n": band_n}

def style_gene_axis(ax):
    ax.set_xlim(XMIN, XMAX)
    ax.set_ylim(0, 4.15)
    for side in ["top", "right", "left"]:
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_visible(True)
    ax.spines["bottom"].set_linewidth(BASELINE_LW)
    ax.spines["bottom"].set_color(AXIS_COLOR)
    ax.set_yticks([])
    xticks = [1, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000]
    ax.set_xticks(xticks)
    ax.tick_params(axis="x", width=TICK_LW, length=2.3,
                   color=AXIS_COLOR, labelsize=6.6,
                   labelcolor=AXIS_TEXT_COLOR, pad=2.5)
    ax.set_xlabel("HXB2 nucleotide position", fontsize=8.2,
                  color="#333333", labelpad=4)

def draw_gene_box(ax, label, start, end, y, h=0.46):
    ax.add_patch(Rectangle(
        (start, y), end - start, h,
        facecolor="white", edgecolor="#444444",
        linewidth=0.48, zorder=2
    ))
    ax.text((start + end) / 2.0, y + h / 2.0, label,
            ha="center", va="center",
            fontsize=GENE_LABEL_FS, color="#333333", zorder=3)

def draw_primary_gene_tracks(ax):
    for label, start, end in PRIMARY_TRACK_1:
        draw_gene_box(ax, label, start, end, 0.47)
    for label, start, end in PRIMARY_TRACK_2:
        draw_gene_box(ax, label, start, end, 1.17)

def draw_split_gene(ax, label, exon1, exon2, y):
    h = 0.40
    s1, e1 = exon1
    s2, e2 = exon2
    ax.add_patch(Rectangle((s1, y), e1 - s1, h,
                           facecolor="white", edgecolor="#444444",
                           linewidth=0.48, zorder=2))
    ax.add_patch(Rectangle((s2, y), e2 - s2, h,
                           facecolor="white", edgecolor="#444444",
                           linewidth=0.48, zorder=2))
    mid_y = y + h / 2.0
    ax.plot([e1, s2], [mid_y, mid_y], color="#555555", lw=0.40, zorder=1)
    ax.text((s1 + e2) / 2.0, y + h + 0.08, label,
            ha="center", va="bottom",
            fontsize=GENE_LABEL_FS, color="#333333")

def draw_split_gene_tracks(ax):
    draw_split_gene(ax, "tat", SPLIT_GENES["tat"][0], SPLIT_GENES["tat"][1], 2.05)
    draw_split_gene(ax, "rev", SPLIT_GENES["rev"][0], SPLIT_GENES["rev"][1], 2.85)

def build_legend_handles():
    handles = []
    for s in SNP_ORDER:
        handles.append(Line2D([0], [0], marker="s", linestyle="None",
                              markerfacecolor=SNP_COLORS[s],
                              markeredgecolor=SNP_COLORS[s],
                              markersize=5.2, label=s))
    handles.extend([
        Line2D([0], [0], marker="D", linestyle="None",
               markerfacecolor="white", markeredgecolor=DEL_EDGE,
               markersize=4.8, label="DEL event"),
        Line2D([0], [0], marker="v", linestyle="None",
               markerfacecolor="white", markeredgecolor=INS_EDGE,
               markersize=5.1, label="INS event"),
        Line2D([0], [0], marker="D", linestyle="None",
               markerfacecolor=NOUNIQUE_COLOR, markeredgecolor=NOUNIQUE_COLOR,
               markersize=4.8, label="No unique consensus"),
        Line2D([0], [0], color=UNRESOLVED_COLOR,
               linewidth=4.0, label="Unresolved region"),
    ])
    return handles

for path in [PLOT_TABLE, META_TABLE, NOUNIQUE_TABLE, PLOT_QC]:
    require(path)

plot_rows = read_tsv(PLOT_TABLE)
meta_rows = read_tsv(META_TABLE)
nq_rows = read_tsv(NOUNIQUE_TABLE)
qc_rows = read_tsv(PLOT_QC)

ready_rows = [r for r in qc_rows if r.get("Check") == "READY_FOR_FIGURE"]
if not ready_rows or ready_rows[0].get("Observed") != "YES":
    raise RuntimeError("Step06 plot-data QC does not report READY_FOR_FIGURE = YES.")

patient_meta = OrderedDict()
for r in sorted(meta_rows, key=lambda x: patient_sort_key(x["Patient"])):
    p = r["Patient"].strip()
    patient_meta[p] = {
        "HXB2_start": safe_int(r.get("HXB2_start"), XMIN),
        "HXB2_end": safe_int(r.get("HXB2_end"), XMAX),
    }

patient_order = list(patient_meta)
if len(patient_order) != 11:
    raise RuntimeError(f"Expected 11 patient panels; found {len(patient_order)}.")

plot_by_patient = defaultdict(list)
for r in plot_rows:
    plot_by_patient[r["Patient"].strip()].append(r)

nq_by_patient = defaultdict(list)
for r in nq_rows:
    nq_by_patient[r["Patient"].strip()].append(r)

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.size": 8,
    "axes.linewidth": 0.4,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
})

n_patients = len(patient_order)
fig = plt.figure(figsize=(15.2, n_patients * 1.16 + 2.45))
gs = fig.add_gridspec(
    nrows=n_patients + 1,
    ncols=1,
    height_ratios=[1.0] * n_patients + [1.65],
    hspace=0.075,
)

max_stacked_height = 0
total_no_unique_points = 0
total_unresolved_bands = 0
event_label_count = 0

for i, patient in enumerate(patient_order):
    ax = fig.add_subplot(gs[i, 0])
    style_patient_axis(ax)

    rows = plot_by_patient.get(patient, [])
    snp_rows = [r for r in rows if r.get("Mutation_type") == "SNP"]
    del_rows = [r for r in rows if r.get("Mutation_type") == "DEL"]
    ins_rows = [r for r in rows if r.get("Mutation_type") == "INS"]

    max_stacked_height = max(max_stacked_height, draw_snp_bars(ax, snp_rows))

    del_placed = draw_event_track(
        ax, del_rows, DEL_Y, "D", DEL_MARKER_SIZE, DEL_LABEL_LEVELS
    )
    ins_placed = draw_event_track(
        ax, ins_rows, INS_Y, "v", INS_MARKER_SIZE, INS_LABEL_LEVELS
    )
    event_label_count += len(del_placed) + len(ins_placed)

    nq_stats = draw_no_unique(ax, nq_by_patient.get(patient, []))
    total_no_unique_points += nq_stats["point_n"]
    total_unresolved_bands += nq_stats["band_n"]

    meta = patient_meta[patient]
    draw_patient_label(ax, patient, meta["HXB2_start"], meta["HXB2_end"])

gene_ax = fig.add_subplot(gs[n_patients, 0])
style_gene_axis(gene_ax)
draw_primary_gene_tracks(gene_ax)
draw_split_gene_tracks(gene_ax)

fig.suptitle(
    "Patient-specific mutation landscape projected to HXB2 coordinates",
    fontsize=TITLE_FS, y=0.992
)

fig.legend(
    handles=build_legend_handles(),
    loc="upper center",
    bbox_to_anchor=(0.55, 0.978),
    ncol=8,
    frameon=False,
    fontsize=5.8,
    handletextpad=0.35,
    columnspacing=0.80,
)

plt.subplots_adjust(left=0.155, right=0.992, top=0.950, bottom=0.055)

fig.savefig(OUT_SVG, bbox_inches="tight")
fig.savefig(OUT_PDF, bbox_inches="tight")
fig.savefig(OUT_PNG, dpi=600, bbox_inches="tight")
plt.close(fig)

qc_out = [
    {"Check": "Patient_panel_n", "Value": len(patient_order),
     "Status": "PASS" if len(patient_order) == 11 else "FAIL"},
    {"Check": "Maximum_stacked_SNP_height", "Value": max_stacked_height, "Status": "INFO"},
    {"Check": "Common_y_axis_max", "Value": YMAX, "Status": "INFO"},
    {"Check": "SNP_bar_width", "Value": BAR_WIDTH, "Status": "INFO"},
    {"Check": "DEL_marker_size", "Value": DEL_MARKER_SIZE, "Status": "INFO"},
    {"Check": "INS_marker_size", "Value": INS_MARKER_SIZE, "Status": "INFO"},
    {"Check": "Event_count_labels_drawn", "Value": event_label_count, "Status": "INFO"},
    {"Check": "No_unique_point_markers_drawn", "Value": total_no_unique_points, "Status": "INFO"},
    {"Check": "Unresolved_hex_bands_drawn", "Value": total_unresolved_bands, "Status": "INFO"},
    {"Check": "SVG_created", "Value": str(OUT_SVG), "Status": "PASS" if OUT_SVG.exists() else "FAIL"},
    {"Check": "PDF_created", "Value": str(OUT_PDF), "Status": "PASS" if OUT_PDF.exists() else "FAIL"},
    {"Check": "PNG_created", "Value": str(OUT_PNG), "Status": "PASS" if OUT_PNG.exists() else "FAIL"},
]

write_tsv(OUT_QC, ["Check", "Value", "Status"], qc_out)

print("=" * 82)
print("FIGURE 1 PATIENT MUTATION LANDSCAPE COMPLETE")
print("=" * 82)
print()
print(f"Patients                         : {len(patient_order)}")
print(f"Maximum stacked SNP height       : {max_stacked_height}")
print(f"Common y-axis maximum            : {YMAX}")
print(f"No-unique point markers          : {total_no_unique_points}")
print(f"Unresolved hexagonal bands       : {total_unresolved_bands}")
print(f"Event count labels               : {event_label_count}")
print()
print("Outputs:")
print(OUT_SVG)
print(OUT_PDF)
print(OUT_PNG)
print(OUT_QC)
