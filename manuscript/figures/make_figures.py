#!/usr/bin/env python3
"""Main figures for the revised manuscript (Genome Research format), from the result tables of the pipeline
(main run: <data>/work/...) and of the revision stages (<rev>/: results/revision tables).

Usage: python make_figures.py <project dir with work/> <results/revision> <out_dir> [figure numbers]
       env RIDL_REV5: results/revision5 (default <project dir>/results/revision5)
Okabe-Ito colours; Liberation Sans (Helvetica/Arial metrics); panel letters 10 pt bold capitals.
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from matplotlib.legend_handler import HandlerTuple
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrow, Patch, Rectangle
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd

D = sys.argv[1] if len(sys.argv) > 1 else "data"
RV = sys.argv[2] if len(sys.argv) > 2 else "rev"
OUT = sys.argv[3] if len(sys.argv) > 3 else "out_v2"
W = os.path.join(D, "work")
R5 = os.environ.get("RIDL_REV5", os.path.join(D, "results", "revision5"))
os.makedirs(OUT, exist_ok=True)

# Liberation Sans ships with the code (manuscript/fonts, SIL OFL 1.1), so figures look the same on any machine
import glob as _glob
from matplotlib import font_manager as _fm
for _f in _glob.glob(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "fonts", "*.ttf")):
    _fm.fontManager.addfont(_f)
plt.rcParams.update({
    "font.family": "Liberation Sans", "font.size": 7.5, "axes.titlesize": 8, "axes.labelsize": 7.5,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 6.8, "axes.linewidth": 0.6,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6, "xtick.major.size": 2.5, "ytick.major.size": 2.5,
    "axes.spines.top": False, "axes.spines.right": False, "pdf.fonttype": 42, "svg.fonttype": "none",
    "axes.titleweight": "bold", "axes.titlelocation": "left", "legend.frameon": False,
})
INK, GRID = "#222222", "#d9d9d9"
# Okabe-Ito
OI = {"orange": "#E69F00", "sky": "#56B4E9", "green": "#009E73", "yellow": "#F0E442", "blue": "#0072B2",
      "vermillion": "#D55E00", "purple": "#CC79A7", "black": "#000000", "grey": "#999999"}
C = {"AluS": OI["orange"], "AluJ": OI["orange"], "AluY": OI["orange"], "MIR": OI["green"], "MIR3": OI["green"],
     "L2abc": OI["purple"], "L2d": OI["purple"], "LTR": OI["sky"], "L1": OI["vermillion"], "ord": "#4D4D4D",
     "pc": OI["blue"], "exon": OI["vermillion"], "intron": "#8C8C8C", "intergenic": "#C9C9C9", "te": OI["vermillion"],
     "acc": OI["vermillion"], "don": OI["blue"]}
LABEL = {"AluS": "AluS", "AluJ": "AluJ", "AluY": "AluY", "MIR": "MIR", "MIR3": "MIR3", "L2abc": "L2a/b/c",
         "L2d": "L2d", "LTR": "LTR", "L1": "L1"}
GROUPS = ["AluS", "AluJ", "AluY", "MIR", "MIR3", "L2abc"]
STAR = {"AluS": "AluS*", "AluJ": "AluJ*", "AluY": "AluY*", "MIR": "MIR*", "MIR3": "MIR3", "L2abc": "L2abc", "L2d": "L2d*"}
MK = {"AluS": "o", "AluJ": "s", "AluY": "^", "MIR": "o", "MIR3": "s", "L2abc": "o", "L2d": "s"}
CREATED = "#CFCFCF"  # neutral grey (not Okabe-Ito): sites created by mutation
MM = 1 / 25.4
PDF_ALL = None


def tsv(*p):
    return pd.read_csv(os.path.join(W, *p), sep="\t")


def rv(name):
    return pd.read_csv(os.path.join(RV, name), sep="\t")


def r5(name):
    return pd.read_csv(os.path.join(R5, name), sep="\t")


def letter(ax, s, x=-0.12, y=1.04):
    ax.text(x, y, s.upper(), transform=ax.transAxes, fontsize=10, fontweight="bold", va="bottom", ha="left",
            color=INK)


def save(fig, name):
    name += os.environ.get("FIGSUFFIX", "")  # e.g. FIGSUFFIX=_v2 gives test renders their own file names
    fig.savefig(os.path.join(OUT, f"{name}.png"), dpi=300, facecolor="white")
    fig.savefig(os.path.join(OUT, f"{name}.pdf"), facecolor="white")
    if PDF_ALL is not None:
        PDF_ALL.savefig(fig, facecolor="white")
    plt.close(fig)


def whole(run, layer="population"):
    w = tsv("prof_v5", run, "whole_element.tsv")
    return w[w.layer == layer].set_index("orientation")


def ladder(ax, rows, xlim, xlabel, ref=0.0, fs=6.5):
    """rows: list of (label, est, lo, hi, color, marker, filled). Off-scale values get an arrow."""
    n = len(rows)
    y = np.arange(n)[::-1]
    ax.axvline(ref, color=GRID, lw=0.8, zorder=0)
    for yi, (lab, e, lo, hi, col, mk, filled) in zip(y, rows):
        if e > xlim[1] or e < xlim[0]:
            right = e > xlim[1]
            edge = xlim[1] if right else xlim[0]
            span = xlim[1] - xlim[0]
            start = edge - 0.30 * span if right else edge + 0.30 * span
            ax.annotate("", xy=(edge, yi), xytext=(start, yi),
                        arrowprops=dict(arrowstyle="-|>", color=col, lw=1.1, mutation_scale=7))
            ax.text(start - 0.02 * span if right else start + 0.02 * span, yi, f"{e:+.2f}".replace("-", "−"),
                    ha="right" if right else "left", va="center", fontsize=fs, color=col)
            continue
        ax.plot([lo, hi], [yi, yi], color=col, lw=1.2, solid_capstyle="butt")
        ax.plot(e, yi, marker=mk, ms=4.2, color=col, mfc=col if filled else "white", mew=1.0, zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels([r[0] for r in rows])
    ax.set_ylim(-0.7, n - 0.3)
    ax.set_xlim(*xlim)
    ax.set_xlabel(xlabel)
    ax.tick_params(axis="y", length=0)


# ----------------------------------------------------------------------------- Figure 1
def fig1():
    fig = plt.figure(figsize=(180 * MM, 95 * MM))
    gs = GridSpec(2, 2, figure=fig, width_ratios=[1.35, 1], height_ratios=[1, 1], wspace=0.42, hspace=0.8,
                  left=0.03, right=0.98, top=0.93, bottom=0.11)
    ax = fig.add_subplot(gs[:, 0])
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    ax.axis("off")
    letter(ax, "a", x=0.0, y=1.0)
    ax.plot([4.5, 95.5], [90, 90], color=C["intron"], lw=1.2, zorder=1, solid_capstyle="butt")
    EXONS = ((4, 20), (42, 60), (82, 96))
    for x0, x1 in EXONS:
        ax.add_patch(Rectangle((x0, 86.0), x1 - x0, 8, facecolor="#FBE3D4", edgecolor=C["exon"], lw=1.0, zorder=2))

    def te(x, y, w, right=True, col=None):
        if col is None:
            inside = any(a <= x and x + w <= b for a, b in EXONS)
            col = C["exon"] if inside else "#6E6E6E"
        if right:
            ax.add_patch(FancyArrow(x, y, w, 0, width=2.6, head_width=4.2, head_length=2.2,
                                    length_includes_head=True, color=col, zorder=3, lw=0))
        else:
            ax.add_patch(FancyArrow(x + w, y, -w, 0, width=2.6, head_width=4.2, head_length=2.2,
                                    length_includes_head=True, color=col, zorder=3, lw=0))
    for args in ((7, 90, 9, False), (45, 90, 9), (85, 90, 8, False), (25, 90, 7), (33, 90, 6, False), (66, 90, 8)):
        te(*args)
    ax.text(89, 80.5, "exon: mature RNA", ha="center", fontsize=6.8, color=C["exon"])
    ax.text(30, 80.5, "intron: spliced out", ha="center", fontsize=6.8, color="#666666")
    ax.text(4, 99, "lncRNA gene with TE copies", fontsize=7.3, fontweight="bold", color=INK)
    ax.text(4, 73, "Copies placed on consensus coordinates", fontsize=7.3, fontweight="bold", color=INK)
    rng = np.random.default_rng(4)
    x0, x1 = 16, 92
    cx = lambda p: x0 + (p - 1) / 299 * (x1 - x0)
    win = (120, 160)
    ax.add_patch(Rectangle((cx(win[0]), 43), cx(win[1]) - cx(win[0]), 27, color="#FCE6CF", zorder=0, lw=0))
    copies = [(True, 10, 300, 67), (True, 60, 250, 62.5), (True, 1, 190, 58),
              (False, 30, 300, 51.5), (False, 1, 260, 47), (False, 90, 300, 42.5)]
    for exonic, a, b, yy in copies:
        col = C["exon"] if exonic else C["intron"]
        ax.plot([cx(a), cx(b)], [yy, yy], color=col, lw=2.6, solid_capstyle="round", zorder=1)
        n = int((b - a) / 300 * 16)
        for p in rng.uniform(a + 4, b - 4, n):
            if exonic and win[0] - 4 <= p <= win[1] + 4:
                continue
            ax.plot(cx(p), yy, "o", ms=1.6, color=INK, zorder=2)
    ax.text(14, 62.5, "exonic", ha="right", va="center", fontsize=6.8, color=C["exon"])
    ax.text(14, 47, "intronic", ha="right", va="center", fontsize=6.8, color="#666666")
    ax.plot([x0, x1], [38.5, 38.5], color=INK, lw=0.6)
    for p in (1, 100, 200, 300):
        ax.plot([cx(p)] * 2, [38.5, 37], color=INK, lw=0.6)
        ax.text(cx(p), 33.5, str(p), ha="center", fontsize=6.3)
    ax.text(x1, 29.5, "consensus position", ha="right", fontsize=6.3, color="#555555")
    ax.plot(cx(252), 73.6, "o", ms=1.6, color=INK)
    ax.text(cx(252) + 1.5, 73.6, "variant", fontsize=6.3, va="center", color="#555555")
    ax.text(4, 25, "Compare, window by window", fontsize=7.3, fontweight="bold", color=INK)
    p = np.arange(1, 301, 5)
    base = 13 + 0.6 * np.sin(p * 0.37)
    dip = np.exp(-0.5 * ((p - 140) / 9) ** 2)
    exo = 13 + 0.6 * np.sin(p * 0.29 + 1) - 7 * dip
    ax.add_patch(Rectangle((cx(win[0]), 3), cx(win[1]) - cx(win[0]), 18, color="#FCE6CF", zorder=0, lw=0))
    ax.plot(cx(p), base, color=C["intron"], lw=1.2)
    ax.plot(cx(p), exo, color=C["exon"], lw=1.4)
    ax.text(cx(win[1]) + 1.5, 6, "fewer variants in exonic copies\n= candidate RNA domain", fontsize=6.3,
            va="center", color=C["exon"])
    ax.text(x0 - 2, 13, "variant\nfrequency", ha="right", va="center", fontsize=6.3, color="#555555")

    axb = fig.add_subplot(gs[0, 1])
    letter(axb, "b", x=-0.30)
    rows = []
    for g in GROUPS:
        w = whole(f"{g}_b20")
        rows.append((g, int(w.loc["sense", "units_fg"]), int(w.loc["antisense", "units_fg"])))
    y = np.arange(len(rows))[::-1]
    for yi, (g, s, a) in zip(y, rows):
        # patch 16: antisense dark, sense pale, as in Fig. 4B
        axb.barh(yi + 0.19, s, height=0.36, color=C[g], alpha=0.45)
        axb.barh(yi - 0.19, a, height=0.36, color=C[g])
        axb.text(max(s, a) + 20, yi, f"{s:,} / {a:,}", va="center", fontsize=6.3)
    axb.set_yticks(y)
    axb.set_yticklabels([LABEL[r[0]] for r in rows])
    axb.tick_params(axis="y", length=0)
    axb.set_xlabel("exonic TE units (copies or their exonic parts)")
    axb.set_xlim(0, 1400)
    axb.legend(handles=[Patch(color="#999999", alpha=0.45, label="sense"), Patch(color="#999999", label="antisense")],
               loc="lower right", ncol=1, handlelength=1, borderaxespad=0)
    axb.set_title("Exonic TE units per group")

    axc = fig.add_subplot(gs[1, 1])
    letter(axc, "c", x=-0.30)
    sizes = {10: "power", 20: "power_b20", 30: "power_b30"}
    for g in GROUPS:
        ys = []
        for s, d in sizes.items():
            t = tsv(d, "power_summary.tsv")
            t = t[t.subfamily == STAR[g]]
            ys.append(100 * t["median_rel_mde"].mean())
        axc.plot(list(sizes), ys, "-", marker=MK[g], color=C[g], ms=3.2, lw=1.1, label=LABEL[g])
    axc.axhline(10, color=INK, lw=0.6, ls="--")
    axc.text(31.0, 9.6, "10%", fontsize=6.3, va="top", ha="left")
    axc.set_xticks([10, 20, 30])
    axc.set_xlim(8, 34)
    axc.set_ylim(0, 31)
    axc.set_xlabel("window along the consensus (nt)")
    axc.set_ylabel("smallest detectable reduction\nin exonic variants (%)")
    axc.legend(ncol=2, loc="upper right", handlelength=1.2, columnspacing=0.8, borderaxespad=0)
    axc.set_title("Smallest detectable reduction (expected)")
    save(fig, "Fig1")


# ----------------------------------------------------------------------------- Figure 2
def fig2():
    fig = plt.figure(figsize=(180 * MM, 150 * MM))
    top = GridSpec(1, 2, figure=fig, width_ratios=[1.5, 0.8], wspace=0.5, left=0.1, right=0.985, top=0.93,
                   bottom=0.6)
    bot = GridSpec(1, 2, figure=fig, width_ratios=[1.75, 1.0], wspace=0.42, left=0.1, right=0.985, top=0.4,
                   bottom=0.115)
    # B: two adjacent contrasts with their own intervals (no gradient implied): intron vs intergenic, exon vs intron
    nig_b = tsv("prof_nonte", "intron_vs_intergenic", "whole_element.tsv").query("layer=='population'").iloc[0]
    nex_b = tsv("prof_nonte", "exon_vs_intron", "whole_element.tsv").query("layer=='population'").iloc[0]
    rows_b = [("repeat-free", C["ord"], True, nig_b, nex_b)]
    for g in ("AluS", "MIR"):
        ig = tsv("diag6", f"{g}_intron_vs_intergenic", "whole_element.tsv").query("layer=='population'")
        ig = ig.set_index("orientation")
        ex = whole(f"{g}_b20")
        for o in ("sense", "antisense"):
            rows_b.append((f"{LABEL[g]} {o}", C[g], o == "sense", ig.loc[o], ex.loc[o]))
    # grouped bars: two separate contrasts per group (no gradient implied); pale = intron vs intergenic, solid = exon vs intron
    ax = fig.add_subplot(top[0, 1])
    w = 0.38
    for k, (lab, col, _, ig, ex) in enumerate(rows_b):
        for d, e, face in ((-w / 2, ig, mcolors.to_rgba(col, 0.35)), (w / 2, ex, col)):
            ax.bar(k + d, e.effect, width=w * 0.95, color=face, edgecolor=col, lw=0.6)
            ax.plot([k + d] * 2, [e.ci_lo, e.ci_hi], color=INK, lw=0.7)
    ax.axhline(0, color=INK, lw=0.6)
    ax.set_xticks(range(len(rows_b)))
    ax.set_xticklabels([r[0].replace("repeat-free", "repeat-\nfree").replace(" sense", "\ns").replace(" antisense", "\nas")
                        for r in rows_b], fontsize=6.6)
    ax.tick_params(axis="x", length=0)
    ax.set_xlim(-0.6, len(rows_b) - 0.4)
    ax.set_ylim(-0.085, 0.012)
    ax.set_ylabel("log OR of a variable site")
    ax.legend(handles=[Patch(facecolor=mcolors.to_rgba("#4D4D4D", 0.35), edgecolor="#4D4D4D", label="intron vs intergenic"),
                       Patch(facecolor="#4D4D4D", edgecolor="#4D4D4D", label="exon vs intron")],
              loc="lower right", fontsize=6.0, borderaxespad=0.1, handlelength=1.2)
    letter(ax, "b", x=-0.36)
    ax.set_title("Introns already carry fewer\nvariants than intergenic DNA")

    axb = fig.add_subplot(top[0, 0])
    letter(axb, "a", x=-0.2)
    rows = []
    for g in GROUPS:
        a, p = whole(f"{g}_b20"), whole(f"{g}_b20_paired")
        for ori, tag in (("sense", "s"), ("antisense", "as")):
            rows.append((f"{LABEL[g]} {tag}", a.loc[ori], p.loc[ori], C[g]))
    o_all = tsv("prof_nonte", "exon_vs_intron", "whole_element.tsv").query("layer=='population'").iloc[0]
    o_par = tsv("prof_nonte", "exon_vs_intron_paired", "whole_element.tsv").query("layer=='population'").iloc[0]
    rows.append(("repeat-free", o_all, o_par, C["ord"]))
    n = len(rows)
    y = np.arange(n)[::-1]
    # patch 9: each row shows the shift from all genes (open) to the same genes (filled); grey band, repeat-free
    # sequence in the same genes. Intervals: 95%, thin for all genes, thick for the same genes.
    axb.axvspan(o_par.ci_lo, o_par.ci_hi, color="#E8E8E8", lw=0, zorder=0)
    axb.axvline(0, color=GRID, lw=0.8, zorder=0)
    for yi, (lab, A, P, col) in zip(y, rows):
        axb.plot([A.ci_lo, A.ci_hi], [yi + 0.2] * 2, color=col, lw=0.7, alpha=0.45)
        axb.plot([P.ci_lo, P.ci_hi], [yi - 0.2] * 2, color=col, lw=1.4)
        axb.annotate("", xy=(P.effect, yi - 0.2), xytext=(A.effect, yi + 0.2),
                     arrowprops=dict(arrowstyle="-|>", color=col, lw=0.8, mutation_scale=6, shrinkA=2.2, shrinkB=2.6))
        axb.plot(A.effect, yi + 0.2, "o", ms=3.4, mfc="white", mec=col, mew=0.9, zorder=3)
        axb.plot(P.effect, yi - 0.2, "o", ms=3.8, color=col, zorder=3)
    axb.set_yticks(y)
    axb.set_yticklabels([r[0] for r in rows], fontsize=6.6)
    axb.tick_params(axis="y", length=0)
    axb.set_xlim(-0.22, 0.28)
    axb.set_ylim(-0.7, n - 0.3)
    axb.set_xlabel("exon vs intron, log OR of a variable site (below 0, fewer variants in exons)")
    axb.legend(handles=[Line2D([], [], marker="o", ls="", mfc="white", mec="#555555", ms=3.4, label="all genes"),
                        Line2D([], [], marker="o", ls="", color="#555555", ms=3.6, label="same genes"),
                        Patch(color="#E8E8E8", label="repeat-free, same genes")],
               loc="upper center", bbox_to_anchor=(0.45, -0.2), ncol=3, fontsize=6.3, handletextpad=0.3,
               columnspacing=0.9, borderaxespad=0)
    axb.set_title("Within genes, most TE groups match repeat-free sequence")

    sub = GridSpecFromSubplotSpec(2, 1, subplot_spec=bot[0, 1], height_ratios=[2.2, 1], hspace=0.12)
    axc = fig.add_subplot(sub[0])
    axd = fig.add_subplot(sub[1], sharex=axc)
    letter(axc, "d", x=-0.32)
    b = tsv("comp", "lnc_paired_bins.tsv")
    x = np.arange(len(b))
    axc.plot(x, b.phylop_p447_exon, "-o", color=C["exon"], ms=3.4, lw=1.2, label="exon")
    axc.plot(x, b.phylop_p447_intron, "-o", color=C["intron"], ms=3.4, lw=1.2, label="intron")
    axc.set_ylabel("mean phyloP")
    axc.legend(loc="upper right", borderaxespad=0)
    plt.setp(axc.get_xticklabels(), visible=False)
    axc.set_title("Lower exonic conservation\nfollows GC content")
    fe = b.units_exon / b.units_exon.sum()
    fi = b.units_intron / b.units_intron.sum()
    axd.bar(x - 0.18, fe, width=0.36, color=C["exon"])
    axd.bar(x + 0.18, fi, width=0.36, color=C["intron"])
    axd.set_ylabel("share of\nwindows", fontsize=6.8)
    axd.set_xticks(x)
    axd.set_xticklabels([r.replace("-", "–") for r in b.gc_range], fontsize=6.2, rotation=35, ha="right")
    axd.set_xlabel("GC content of 100-nt windows")

    # C: decomposition of the all-genes deficit
    axe = fig.add_subplot(bot[0, 0])
    letter(axe, "c", x=-0.2, y=1.2)
    dte, dnt = rv("dropout_te_decomposition.tsv"), rv("dropout_nonte_decomposition.tsv")
    rows = []
    for g in ("AluS", "AluJ", "MIR", "L2abc"):
        for ori, tag in (("sense", "s"), ("antisense", "as")):
            t = dte[(dte.group == STAR[g]) & (dte.orientation == ori)].set_index("design")
            rows.append((f"{LABEL[g]} {tag}", t, C[g]))
    rows.append(("repeat-free", dnt.set_index("design"), C["ord"]))
    # patch 9: grouped bars, as in panel B: solid, exonic copies (genes without intronic copies vs genes with both);
    # pale, intronic copies (genes without exonic copies vs genes with both); lines, 95% intervals
    w = 0.38
    axe.axhline(0, color=INK, lw=0.6)
    for k, (lab, t, col) in enumerate(rows):
        for d, des, face in ((-w / 2, "exon unpaired vs exon paired", col),
                             (w / 2, "intron unpaired vs intron paired", mcolors.to_rgba(col, 0.35))):
            r = t.loc[des]
            axe.bar(k + d, r.mh_log_or, width=w * 0.95, color=face, edgecolor=col, lw=0.6)
            axe.plot([k + d] * 2, [r.ci_lo, r.ci_hi], color=INK, lw=0.7)
    axe.set_xticks(range(len(rows)))
    axe.set_xticklabels([r[0].replace("repeat-free", "repeat-\nfree").replace(" s", "\ns").replace(" as", "\nas")
                         if r[0] != "repeat-free" else "repeat-\nfree" for r in rows], fontsize=6.6)
    axe.tick_params(axis="x", length=0)
    axe.set_xlim(-0.6, len(rows) - 0.4)
    axe.set_ylabel("log OR of a variable site,\ngenes dropping out vs genes kept")
    axe.legend(handles=[Patch(facecolor="#4D4D4D", edgecolor="#4D4D4D",
                              label="exonic copies: genes without intronic copies vs genes with both"),
                        Patch(facecolor=mcolors.to_rgba("#4D4D4D", 0.35), edgecolor="#4D4D4D",
                              label="intronic copies: genes without exonic copies vs genes with both")],
               loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=1, fontsize=6.4, borderaxespad=0.3,
               handletextpad=0.4, handlelength=1.2)
    axe.set_title("The deficit across genes comes mostly from the genes that drop out", pad=30)
    save(fig, "Fig2")


# ----------------------------------------------------------------------------- Figure 3
def fig3():
    fig = plt.figure(figsize=(180 * MM, 150 * MM))
    gs = GridSpec(2, 2, figure=fig, width_ratios=[1, 1], height_ratios=[1, 1.05], wspace=0.5, hspace=0.62,
                  left=0.17, right=0.92, top=0.925, bottom=0.08)
    # A: GC-matched exon - intron conservation in the same genes, mammals (left) and primates (right)
    subA = GridSpecFromSubplotSpec(1, 2, subplot_spec=gs[0, 0], wspace=0.2)
    axA = []

    def gm(path, track):
        t = pd.read_csv(path, sep="\t").set_index(["metric", "comparison"]).loc[(track, "gc_matched")]
        return t.exon_minus_intron, t.ci_lo, t.ci_hi
    srcA = [("protein-coding", os.path.join(W, "comp", "pc_paired_summary.tsv"), C["pc"], "s"),
            ("lncRNA, repeat-free", os.path.join(W, "comp", "lnc_paired_summary.tsv"), C["ord"], "s")]
    for g in ("AluS", "AluJ", "MIR", "L2abc"):
        srcA.append((f"lncRNA, {LABEL[g]}", os.path.join(W, "comp", f"te_{g}_paired_summary.tsv"), C[g], "o"))
    srcA.append(("lncRNA, LTR", os.path.join(RV, "ltr_comp_te_LTR_paired_summary.tsv"), C["LTR"], "o"))
    srcA.append(("lncRNA, L1", os.path.join(RV, "l1lin_comp_te_L1_paired_summary.tsv"), C["L1"], "o"))
    for k, (track, sub_title) in enumerate((("phylop_p447", "447 mammals"), ("phylop_prim", "primates"))):
        ax = fig.add_subplot(subA[0, k])
        if k == 0:
            letter(ax, "a", x=-1.08)
        rows = [(lab, *gm(path, track), col, mk, True) for lab, path, col, mk in srcA]
        ladder(ax, rows, (-0.1, 0.16), "")
        ax.set_xticks([-0.1, 0.0, 0.1])
        if k:
            ax.set_yticklabels([])
        ax.set_title(sub_title, fontsize=7, fontweight="normal", loc="center")
        axA.append(ax)
    p0, p1 = axA[0].get_position(), axA[1].get_position()
    fig.text(p0.x0, 0.993, "GC-matched conservation, same genes", fontsize=8, fontweight="bold", ha="left", va="top")
    fig.text((p0.x0 + p1.x1) / 2, p0.y0 - 0.043, "exon − intron, mean phyloP", fontsize=7.5, ha="center", va="top")

    axb = fig.add_subplot(gs[0, 1])
    letter(axb, "b", x=-0.36)
    B = rv("bounds_bounds.tsv")
    B = B[(B.table == "main") & (B.layer == "excess_poly_noncpg")]
    LB = {"LTR": rv("ltr_bounds_class_bounds.tsv"), "L1": rv("l1lin_bounds_class_bounds.tsv")}
    order = ["AluS", "AluJ", "MIR", "L2abc", "LTR", "L1"]

    def get(g, ori):
        t = LB[g] if g in LB else B[B.subfamily == STAR[g]]
        t = t[(t.layer == "excess_poly_noncpg") & (t.orientation == ori)].iloc[0]
        return t.mean_pct, t.ci_lo_pct, t.ci_hi_pct
    # patch 8: neutral grey "excluded" band, limits as black numbers in a right-hand column,
    # orientations in grey, colour only for the pooled estimate
    y = np.arange(len(order))[::-1] * 1.0
    axb.axvline(0, color="#9a9a9a", lw=0.8, zorder=0)
    for yi, g in zip(y, order):
        e, lo, hi = get(g, "pooled")
        axb.add_patch(Rectangle((-14, yi - 0.36), lo + 14, 0.72, color="#ececec", lw=0, zorder=0))
        for off, ori, mk in ((0.22, "sense", "^"), (-0.22, "antisense", "v")):
            e2, lo2, hi2 = get(g, ori)
            axb.plot([lo2, hi2], [yi + off] * 2, color="#a6a6a6", lw=0.7, zorder=1)
            axb.plot(e2, yi + off, mk, color="#8c8c8c", ms=2.8, zorder=2)
        axb.plot([lo, hi], [yi, yi], color=C[g], lw=1.8, zorder=3)
        axb.plot(e, yi, "o", color=C[g], ms=4.8, zorder=4, mec="white", mew=0.5)
        axb.text(10.6, yi, f"{abs(lo):.1f}%", ha="left", va="center", fontsize=7, color=INK)
    axb.text(10.6, y[0] + 0.55, "excluded", ha="left", va="bottom", fontsize=6.3, color="#555555")
    axb.set_yticks(y)
    axb.set_yticklabels([LABEL[g] for g in order])
    axb.tick_params(axis="y", length=0)
    axb.set_xlim(-14, 10)
    axb.set_ylim(-1.05, len(order) + 0.75)
    axb.set_xlabel("TE-specific change in variant frequency (%)")
    axb.text(-13.5, -0.98, "fewer variants \u2190", fontsize=6.3, color="#555555", va="bottom")
    axb.text(9.5, -0.98, "\u2192 more variants", fontsize=6.3, color="#555555", va="bottom", ha="right")
    axb.legend(handles=[Line2D([], [], marker="o", ls="-", color="#555555", ms=4, label="both orientations"),
                        Line2D([], [], marker="^", ls="-", lw=0.7, color="#9a9a9a", ms=3, label="sense"),
                        Line2D([], [], marker="v", ls="-", lw=0.7, color="#9a9a9a", ms=3, label="antisense"),
                        Patch(color="#ececec", label="excluded at 95%")],
               loc="upper left", ncol=2, fontsize=6.2, handlelength=1.3,
               columnspacing=0.8, handletextpad=0.3, borderaxespad=0.3,
               frameon=True, facecolor="white", edgecolor="none", framealpha=1.0)
    axb.set_title("Limits on TE-specific constraint")

    axc = fig.add_subplot(gs[1, 0])
    letter(axc, "c", x=-0.55, y=1.12)
    pa, po = rv("posctrl_alu_prim_q0.99_posctrl.tsv"), rv("posctrl_old_p447_t1.0_posctrl.tsv")
    rows = []
    for T, g, note in ((pa, "AluS", "primates"), (pa, "AluJ", "primates"), (pa, "AluY", "primates"),
                       (po, "MIR", "mammals"), (po, "L2abc", "mammals"), (po, "L2d", "mammals"), (po, "MIR3", "mammals")):
        t = T[(T.group == STAR[g]) & (T.compartment == "all")]
        A = t[t.comparison.str.startswith("A")].iloc[0]
        Bw = t[t.comparison.str.startswith("B")].iloc[0]
        rows.append((f"{LABEL[g]} ({note})", A, Bw, C[g]))
    y = np.arange(len(rows))[::-1]
    pct = lambda v: 100 * (1 - np.exp(v))
    axc.axvline(0, color=GRID, lw=0.8, zorder=0)
    for yi, (lab, A, Bw, col) in zip(y, rows):
        axc.barh(yi + 0.19, pct(A.mh_log_or), height=0.36, color=col, lw=0)
        axc.plot([pct(A.ci_hi), pct(A.ci_lo)], [yi + 0.19] * 2, color="#222222", lw=0.9)
        axc.barh(yi - 0.19, pct(Bw.mh_log_or), height=0.36, color=col, alpha=0.4, lw=0)
        axc.plot([pct(Bw.ci_hi), pct(Bw.ci_lo)], [yi - 0.19] * 2, color="#222222", lw=0.9)
    axc.set_yticks(y)
    axc.set_yticklabels([r[0] for r in rows])
    axc.tick_params(axis="y", length=0)
    axc.set_ylim(-0.7, len(rows) - 0.3)
    axc.set_xlim(-8, 16)
    axc.set_xlabel("fewer human variants in conserved TE stretches (%)")
    axc.legend(handles=[Patch(color="#555555", label="vs same positions, other copies"),
                        Patch(color="#555555", alpha=0.4, label="vs rest of the same copies")],
               loc="lower left", bbox_to_anchor=(0.0, 1.0), ncol=2, fontsize=6.2, borderaxespad=0.2,
               handletextpad=0.3, columnspacing=1.0)

    # D: planted constraint, AluS and MIR: domain in half of the exonic copies (full grid: Supplementary Fig. S3)
    axc.set_title("Positive control: conserved TE stretches", pad=16)
    sub = GridSpecFromSubplotSpec(2, 2, subplot_spec=gs[1, 1], hspace=0.6, wspace=0.12)
    P = pd.concat([rv(f"power_{g}_{o}_power.tsv") for g in ("AluS", "MIR") for o in ("sense", "antisense")])
    reds = (0.1, 0.25, 0.5)
    im = None
    for i, g in enumerate(("AluS*", "MIR*")):
        for j, o in enumerate(("sense", "antisense")):
            axh = fig.add_subplot(sub[i, j])
            if i == 0 and j == 0:
                letter(axh, "d", x=-0.75, y=1.32)
            t = P[(P.group == g) & (P.orientation == o) & (P.reduction > 0) & (P.fraction_units == 0.5)]
            segs = sorted(t.segment.unique(), key=lambda s_: (s_ == "all", int(s_.split("-")[1]) - int(s_.split("-")[0]) if s_ != "all" else 0))
            M = np.full((len(segs), len(reds)), np.nan)
            for si, sg in enumerate(segs):
                for ci, r in enumerate(reds):
                    q = t[(t.segment == sg) & (t.reduction == r)]
                    if len(q):
                        M[si, ci] = q.power_whole.iloc[0] if sg == "all" else q.power_window.iloc[0]
            im = axh.imshow(M, cmap="cividis", vmin=0, vmax=1, aspect="auto")
            for si in range(M.shape[0]):
                for ci in range(M.shape[1]):
                    v = M[si, ci]
                    if np.isfinite(v):
                        axh.text(ci, si, f"{v:.2f}", ha="center",
                                 va="center", fontsize=6.2, color="black" if v > 0.55 else "white")
            nts = [f"{int(sg.split('-')[1]) - int(sg.split('-')[0]) + 1} nt" if sg != "all" else "whole" for sg in segs]
            axh.set_yticks(range(len(segs)))
            axh.set_yticklabels(nts if j == 0 else [""] * len(segs), fontsize=6.5)
            axh.set_xticks(range(len(reds)))
            axh.set_xticklabels([f"−{int(100 * r)}%" for r in reds], fontsize=6.3)
            for sp in axh.spines.values():
                sp.set_visible(False)
            axh.tick_params(length=0)
            axh.set_title(f"{LABEL[g.rstrip('*')]} {o}", fontsize=7, fontweight="normal")
            if i == 1 and j == 0:
                axh.set_xlabel("variants removed in the domain", fontsize=6.6, x=1.06)
    cax = fig.add_axes([0.94, 0.10, 0.008, 0.25])
    cb = fig.colorbar(im, cax=cax)
    cb.ax.tick_params(labelsize=6, length=1.5)
    cb.outline.set_linewidth(0.4)
    cb.set_label("detection rate", fontsize=6.5)
    fig.text(0.585, 0.475, "Detection of planted constraint", fontsize=8, fontweight="bold", ha="left")
    fig.text(0.585, 0.455, "domain in half of the exonic copies", fontsize=6.6, ha="left", color="#444444")
    save(fig, "Fig3")


# ----------------------------------------------------------------------------- Figure 4
def fig4():
    fig = plt.figure(figsize=(180 * MM, 120 * MM))
    gs = GridSpec(2, 2, figure=fig, width_ratios=[1.25, 1], height_ratios=[1, 1], wspace=0.35, hspace=0.62,
                  left=0.07, right=0.98, top=0.9, bottom=0.08)
    ax = fig.add_subplot(gs[0, 0])
    letter(ax, "a", x=-0.12)
    af = tsv("side", "architecture_features.tsv").set_index(["gene_set", "feature"])
    global LOC
    LOC = r5("feature_local_local.tsv")
    feats = [("bg_exon", "exon\nsequence"), ("TSS", "start"), ("donor", "donor"), ("acceptor", "acceptor"), ("TES", "end")]
    classes = [("Alu", OI["orange"]), ("MIR", OI["green"]), ("L2", OI["purple"]), ("L1", OI["vermillion"]),
               ("LTR", OI["sky"]), ("DNA", OI["yellow"])]
    x = np.arange(len(feats))
    for gset, off in (("lnc", -0.2), ("pc", 0.2)):
        bottom = np.zeros(len(feats))
        for cl, col in classes:
            v = np.array([af.loc[(gset, f), cl] for f, _ in feats]) * 100
            ax.bar(x + off, v, width=0.36, bottom=bottom, color=col, edgecolor="white", lw=0.3,
                   alpha=1.0 if gset == "lnc" else 0.5)
            bottom += v
        other = np.array([af.loc[(gset, f), "any TE"] for f, _ in feats]) * 100 - bottom
        ax.bar(x + off, other, width=0.36, bottom=bottom, color="#BBBBBB", edgecolor="white", lw=0.3,
               alpha=1.0 if gset == "lnc" else 0.5)
        # patch 17: TE share of exonic positions 100-500 nt from the feature (step 38); positions 20-100 nt away
        # mostly lie in the same TE copy as the feature, so they are not an independent expectation
        for xi, (f, _), t in zip(x + off, feats, bottom + other):
            top = t
            if f != "bg_exon":
                z = LOC[(LOC.gene_set == gset) & (LOC.feature == f) & (LOC.side == "exon side")
                        & (LOC.falls_in == "exon") & LOC.dist_bin.isin(["100-200", "200-500"])]
                v = 100 * (z["any TE"] * z.n).sum() / z.n.sum()
                ax.plot([xi - 0.2, xi + 0.2], [v, v], color=INK, lw=1.3, solid_capstyle="butt", zorder=5)
                top = max(t, v)
            ax.text(xi, top + 0.9, f"{t:.0f}", ha="center", fontsize=6.2)
    ax.set_xticks(x)
    ax.set_xticklabels([f[1] for f in feats], fontsize=6.8)
    ax.set_ylabel("inside a TE (%)")
    ax.set_ylim(0, 50)
    ax.set_xlabel("left bar: lncRNA · right (pale): protein-coding", fontsize=6.6, color="#555555")
    ax.legend(handles=[Patch(color=c, label=l) for l, c in classes] + [Patch(color="#BBBBBB", label="other")],
              loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=7, fontsize=6.3, handlelength=0.9,
              columnspacing=0.6, handletextpad=0.3, borderaxespad=0.2)
    ax.plot([0.62, 0.9], [48.3, 48.3], color=INK, lw=1.3, solid_capstyle="butt")
    ax.text(0.98, 48.3, "expected: TE share of exonic sequence 100–500 nt away", ha="left", va="center",
            fontsize=6.2, color=INK)
    ax.set_title("TE share of lncRNA starts, splice sites and ends", pad=18)

    axc = fig.add_subplot(gs[0, 1])
    letter(axc, "b", x=-0.22)
    U = r5("boundary_null_summary.tsv")
    gl = ["AluS", "AluJ", "AluY", "MIR", "MIR3", "L2abc"]
    y = np.arange(len(gl))[::-1]
    for yi, g in zip(y, gl):
        for ori, off, alpha in (("antisense", 0.19, 1.0), ("sense", -0.19, 0.45)):
            r = U[(U.group == g) & (U.orientation == ori)].iloc[0]
            v = 100 * r.observed
            axc.barh(yi + off, v, height=0.36, color=C[g], alpha=alpha)
            e = 100 * r.expected
            axc.plot([e, e], [yi + off - 0.17, yi + off + 0.17], color=INK, lw=1.4, solid_capstyle="butt")
            axc.text(max(v, e) + 1.2, yi + off, f"{v:.0f}%", va="center", fontsize=6.2)
    axc.set_yticks(y)
    axc.set_yticklabels([LABEL[g] for g in gl])
    axc.tick_params(axis="y", length=0)
    axc.set_xlim(0, 100)
    axc.set_xlabel("exonic copies containing a splice site of the gene (%)")
    axc.legend(handles=[Patch(color="#999999", label="antisense"), Patch(color="#999999", alpha=0.45, label="sense"),
                        Line2D([], [], color=INK, lw=1.4, marker="|", ms=6, ls="", label="expected by chance")],
               loc="lower right", borderaxespad=0, handlelength=1)
    axc.set_title("Antisense Alu copies often contain a splice site")

    sub = GridSpecFromSubplotSpec(1, 2, subplot_spec=gs[1, :], wspace=0.18)
    hs = tsv("side", "architecture_hotspots.tsv")
    hs = hs[(hs.gene_set == "lnc") & (hs.rfamily == "Alu")]
    # positions used by Alu-derived exons of protein-coding genes, antisense Alu only (Sorek et al. 2002;
    # Lev-Maor et al. 2003); 275 and 279 are drawn as one cluster
    lit = {"antisense": {"acceptor": [116, 255, 277], "donor": [23, 158]},
           "sense": {"acceptor": [], "donor": []}}
    for ci, ori in enumerate(("antisense", "sense")):
        axh = fig.add_subplot(sub[0, ci])
        if ci == 0:
            letter(axh, "c", x=-0.1, y=1.08)
        bins = np.arange(0, 316, 3)
        for feat, sign, col in (("acceptor", 1, C["acc"]), ("donor", -1, C["don"])):
            d = hs[(hs.orientation == ori) & (hs.feature == feat)]
            h, _ = np.histogram(d.cons_pos, bins=bins, weights=d.n)
            axh.bar(bins[:-1] + 1.5, sign * h, width=3, color=col, lw=0)
            top = (h.max() if len(h) else 1) * 1.08
            for p in lit[ori][feat]:
                axh.plot(p, sign * top, marker="v" if sign > 0 else "^", ms=3.8, color=INK)
        axh.axhline(0, color=INK, lw=0.5)
        axh.set_xlim(0, 315)
        axh.set_xlabel("Alu consensus position")
        yl = axh.get_ylim()
        axh.set_yticks(axh.get_yticks())
        axh.set_yticklabels([f"{abs(int(t))}" for t in axh.get_yticks()])
        axh.set_ylim(yl)
        if ci == 0:
            axh.set_ylabel("lncRNA splice sites")
        axh.set_title(f"Splice sites in Alu copies {ori} to the transcript", fontsize=7.5)
        if ci == 1:
            axh.legend(handles=[Patch(color=C["acc"], label="acceptors"), Patch(color=C["don"], label="donors"),
                                (Line2D([], [], marker="v", ls="", color=INK, ms=3.8),
                                 Line2D([], [], marker="^", ls="", color=INK, ms=3.8))],
                       labels=["acceptors", "donors", "sites of Alu exons in protein-\ncoding genes (antisense Alu)"],
                       handler_map={tuple: HandlerTuple(ndivide=None, pad=0.2)},
                       loc="lower right", fontsize=6.3, handlelength=1.6, borderaxespad=0.2)
    save(fig, "Fig4")


# ----------------------------------------------------------------------------- Figure 5
SPECIES = [("panTro6", "chimpanzee"), ("gorGor6", "gorilla"), ("ponAbe3", "orangutan"), ("rheMac10", "macaque"),
           ("calJac4", "marmoset"), ("mm39", "mouse")]


def fig5():
    A = rv("splice_age_summary.tsv")
    g = lambda s, m, sub="all": A[(A.set == s) & (A.measure == m) & (A.subset == sub)].iloc[0]
    base = None
    bpath = os.path.join(RV, "splice_age_base_baseline.tsv")
    if os.path.exists(bpath):
        base = pd.read_csv(bpath, sep="\t")
    fig = plt.figure(figsize=(180 * MM, 135 * MM))
    gs = GridSpec(2, 2, figure=fig, width_ratios=[1, 1.15], height_ratios=[1.25, 1], wspace=0.42, hspace=0.5,
                  left=0.15, right=0.98, top=0.93, bottom=0.09)
    # A: preformed share by subset
    ax = fig.add_subplot(gs[0, 0])
    letter(ax, "a", x=-0.42, y=1.1)
    subs = [("all", "all TE-derived"), ("site_type=acceptor", "acceptors"), ("site_type=donor", "donors"),
            ("rep_family=Alu", "Alu"), ("rep_family=MIR", "MIR"), ("rep_family=L2", "L2"),
            ("te_orientation=antisense", "antisense"), ("te_orientation=sense", "sense")]
    y = np.arange(len(subs))[::-1]
    for yi, (s, lab) in zip(y, subs):
        r = g("lnc_TE_site", "preformed", s)
        ax.barh(yi, 100 * r.share, color=C["te"], height=0.6)
        ax.barh(yi, 100 * (1 - r.share), left=100 * r.share, color=CREATED, height=0.6)
        ax.text(2, yi, f"{100 * r.share:.0f}%", va="center", fontsize=6.3, color="white", fontweight="bold")
    ax.set_yticks(y)
    ax.set_yticklabels([l for _, l in subs])
    ax.tick_params(axis="y", length=0)
    ax.set_xlim(0, 100)
    ax.set_xlabel("TE-derived splice sites (%)")
    ax.legend(handles=[Patch(color=C["te"], label="preformed in the consensus"),
                       Patch(color=CREATED, label="created by mutation")],
              loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=2, fontsize=6.3, borderaxespad=0.2, handlelength=1)
    ax.set_title("About half of TE-derived sites were preformed", pad=16)

    # B: presence by species, splice sites vs function-free dinucleotides (neutral baseline)
    gb = GridSpecFromSubplotSpec(2, 1, subplot_spec=gs[0, 1], height_ratios=[1.55, 1], hspace=0.12)
    axb = fig.add_subplot(gb[0])
    axd = fig.add_subplot(gb[1], sharex=axb)
    letter(axb, "b", x=-0.2, y=1.19)
    xs = np.arange(len(SPECIES))
    for st, lab, nlab, col, dx in (("pc_site", "protein-coding sites", "", C["pc"], -0.12),
                                   ("lnc_nonTE_site", "lncRNA sites outside repeats", "", C["ord"], 0.0),
                                   ("lnc_TE_site", "TE-derived lncRNA sites", "", C["te"], 0.12)):
        t = base[(base.set == st) & (base.condition == "all aligned") & (base.subset == "all")].set_index("species")
        vs = [100 * t.loc[sp, "site_present"] for sp, _ in SPECIES]
        vc = [100 * t.loc[sp, "ctrl_present"] for sp, _ in SPECIES]
        # sites solid, controls dotted and open; series shifted a little sideways so overlapping points stay visible
        axb.plot(xs + dx, vs, "-", marker="o", ms=3.0, color=col, lw=1.3, zorder=3)
        axb.plot(xs + dx, vc, ":", marker="^", ms=3.0, mfc="white", mew=0.8, color=col, lw=1.0, zorder=2)
        if st != "pc_site":
            dd = [100 * t.loc[sp, "diff"] for sp, _ in SPECIES]
            lo = [100 * t.loc[sp, "diff_lo"] for sp, _ in SPECIES]
            hi = [100 * t.loc[sp, "diff_hi"] for sp, _ in SPECIES]
            sh = dx
            axd.errorbar(xs + sh, dd, yerr=[np.subtract(dd, lo), np.subtract(hi, dd)], fmt="o", color=col, ms=3.0,
                         elinewidth=0.8, capsize=0)
        else:
            pcd = [100 * t.loc[sp, "diff"] for sp, _ in SPECIES]
    axd.axhline(0, color=GRID, lw=0.8, zorder=0)
    axd.set_ylim(-4, 12.5)
    axd.set_ylabel("site − control\n(points)")
    axd.text(0.02, 0.95, f"protein-coding: +{pcd[0]:.1f} to +{pcd[-1]:.0f} (off scale)", transform=axd.transAxes,
             fontsize=6.0, va="top", color=C["pc"])
    axd.set_xticks(xs)
    axd.set_xticklabels([n for _, n in SPECIES], rotation=30, ha="right")
    plt.setp(axb.get_xticklabels(), visible=False)
    axb.set_ylabel("present at the orthologous\nposition (% of aligned)")
    axb.set_ylim(40, 101)
    hs = [Line2D([], [], color=C["pc"], lw=1.3, marker="o", ms=3.0, label="protein-coding"),
          Line2D([], [], color=C["ord"], lw=1.3, marker="o", ms=3.0, label="lncRNA, outside repeats"),
          Line2D([], [], color=C["te"], lw=1.3, marker="o", ms=3.0, label="lncRNA, TE-derived"),
          Line2D([], [], color=INK, lw=1.3, marker="o", ms=3.0, label="splice sites"),
          Line2D([], [], color=INK, lw=1.0, ls=":", marker="^", ms=3.0, mfc="white", mew=0.8,
                 label="controls (unused GT/AG)")]
    axb.legend(handles=hs, loc="lower left", fontsize=5.8, borderaxespad=0, handlelength=2.0, labelspacing=0.35)
    axb.set_title("lncRNA splice sites are kept about as\noften as unused dinucleotides")

    # C: age classes
    axc = fig.add_subplot(gs[1, :])
    letter(axc, "c", x=-0.17, y=1.1)
    classes = [(("age=human-specific", "age=shared with panTro6", "age=shared with gorGor6",
                 "age=shared with ponAbe3"), "human and great apes only"),
               (("age=shared with rheMac10",), "up to macaque"), (("age=shared with calJac4",), "up to marmoset"),
               (("age=shared with mm39",), "up to mouse")]
    cc = [OI["vermillion"], OI["orange"], OI["sky"], OI["blue"]]
    tcol = ["white", INK, INK, "white"]
    sets = [("lnc_TE_site", "all", "TE-derived (all)"), ("lnc_TE_site", "rep_family=Alu", "· Alu"),
            ("lnc_TE_site", "rep_family=MIR", "· MIR"), ("lnc_TE_site", "rep_family=L2", "· L2"),
            ("lnc_TE_site", "ancestral=created", "· created"), ("lnc_TE_site", "ancestral=preformed", "· preformed"),
            ("lnc_nonTE_site", "all", "lncRNA, outside\nrepeats"), ("pc_site", "all", "protein-coding")]
    y = np.arange(len(sets))[::-1]
    for yi, (s, sub, lab) in zip(y, sets):
        left = 0
        for (ms, _), col, tc in zip(classes, cc, tcol):
            r = A[(A.set == s) & A.measure.isin(ms) & (A.subset == sub)]
            v = 100 * r.share.sum() if len(r) else 0
            axc.barh(yi, v, left=left, color=col, height=0.66, edgecolor="white", lw=0.6)
            if v >= 4:
                axc.text(left + v / 2, yi, f"{v:.0f}", ha="center", va="center", fontsize=6.2, color=tc)
            left += v
    axc.set_yticks(y)
    axc.set_yticklabels([l for _, _, l in sets])
    axc.tick_params(axis="y", length=0)
    axc.set_xlim(0, 100)
    axc.set_xlabel("splice sites (%), by the most distant species sharing the site")
    axc.legend(handles=[Patch(color=c, label=l) for (_, l), c in zip(classes, cc)], loc="lower center",
               bbox_to_anchor=(0.5, 1.0), ncol=4, fontsize=6.5, borderaxespad=0.2, handlelength=1, columnspacing=1.2)
    axc.set_title("TE-derived sites, especially in Alu, are the youngest class", pad=16)
    save(fig, "Fig6")


# ----------------------------------------------------------------------------- Figure 6
def fig6():
    fig = plt.figure(figsize=(180 * MM, 125 * MM))
    gs = GridSpec(1, 3, figure=fig, wspace=0.22, left=0.25, right=0.985, top=0.92, bottom=0.58)
    st = tsv("splice", "res_stats.tsv")
    si = rv("splice_intron_res_stats.tsv")
    te = "lncRNA TE-derived splice sites vs same positions in unused copies"
    pick = [(st, "protein-coding splice sites vs intronic flank", "all", "protein-coding (vs intron)", C["pc"]),
            (st, "lncRNA non-TE splice sites vs intronic flank", "all", "lncRNA, outside repeats (vs intron)", C["ord"]),
            (st, te, "all", "lncRNA, in TEs (vs unused copies)", C["te"]),
            (si, te, "all", "lncRNA, in TEs (vs intronic copies)", C["te"]),
            (st, te, "Alu", "· Alu", C["AluS"]), (st, te, "MIR", "· MIR", C["MIR"]), (st, te, "L2", "· L2", C["L2abc"])]
    rows_v, rows_s, rows_p = [], [], []
    for T, comp, sub, lab, col in pick:
        r = T[(T.comparison == comp) & (T.subset == sub)].iloc[0]
        mk = "s" if "intronic copies" in lab else "o"
        rows_v.append((lab, r.mh_log_or, r.mh_log_or_lo, r.mh_log_or_hi, col, mk, True))
        rows_s.append((lab, r.singleton_diff, r.singleton_diff_lo, r.singleton_diff_hi, col, mk, True))
        rows_p.append((lab, r.phylop_p447_diff, r.phylop_p447_diff_lo, r.phylop_p447_diff_hi, col, mk, True))
    for i, (rows, xl, xlab, title) in enumerate((
            (rows_v, (-0.5, 0.3), "log OR of a variable site\n(below 0 = fewer variants)", "Human variants"),
            (rows_s, (-0.1, 0.15), "difference in singleton share\n(above 0 = rarer variants)", "Rare variants"),
            (rows_p, (-0.25, 0.45), "difference in phyloP\n(above 0 = more conserved)", "Conservation"))):
        ax = fig.add_subplot(gs[0, i])
        letter(ax, "a" if i == 0 else "", x=-1.1 if i == 0 else 0, y=1.04)
        ladder(ax, rows, xl, xlab)
        if i:
            ax.set_yticklabels([])
        ax.set_title(title)

    # B: by usage tertile (variants and phyloP); C: created vs preformed
    S = rv("splice_strata_total_reads_strata.tsv")
    mp = {'lncRNA non-TE splice sites vs intronic flank': "nonTE", 'protein-coding splice sites vs intronic flank': "pc",
          'lncRNA TE-derived splice sites vs unused copies (main)': "TE"}
    S["c"] = S.comparison.map(mp)
    sub = GridSpec(1, 3, figure=fig, left=0.09, right=0.985, top=0.37, bottom=0.1, wspace=0.5,
                   width_ratios=[1, 1, 0.9])
    sets = [("pc", "protein-coding", C["pc"]), ("nonTE", "lncRNA, outside repeats", C["ord"]),
            ("TE", "lncRNA, in TEs", C["te"])]
    for k, (col_e, lo_, hi_, ylab, title) in enumerate((
            ("mh_log_or", "mh_log_or_lo", "mh_log_or_hi", "log OR of a variable site", "Variants by usage"),
            ("phylop_p447_diff", "phylop_p447_diff_lo", "phylop_p447_diff_hi", "difference in phyloP",
             "Conservation by usage"))):
        ax = fig.add_subplot(sub[0, k])
        if k == 0:
            letter(ax, "b", x=-0.3, y=1.06)
        for ci, (c, lab, col) in enumerate(sets):
            if k == 1 and c == "pc":
                v = [S[(S.c == c) & (S.subset == f"usage q{i}")].iloc[0][col_e] for i in (1, 2, 3)]
                ax.text(0.03, 0.97, f"protein-coding (off scale):\n+{v[0]:.1f} → +{v[1]:.1f} → +{v[2]:.1f}",
                        transform=ax.transAxes, fontsize=6.3, va="top", color=C["pc"])
                continue
            v = [S[(S.c == c) & (S.subset == f"usage q{i}")].iloc[0] for i in (1, 2, 3)]
            xs = np.array([1, 2, 3]) + (ci - 1) * 0.12
            ax.errorbar(xs, [r[col_e] for r in v], yerr=[[r[col_e] - r[lo_] for r in v], [r[hi_] - r[col_e] for r in v]],
                        fmt="-o", color=col, ms=3.4, lw=1.1, elinewidth=0.8, capsize=0, label=lab)
        ax.axhline(0, color=GRID, lw=0.8, zorder=0)
        ax.set_xticks([1, 2, 3])
        ax.set_xticklabels(["low", "middle", "high"])
        ax.set_xlim(0.6, 3.4)
        ax.set_xlabel("junction usage in GTEx (tertile)")
        ax.set_ylabel(ylab)
        ax.set_title(title)
        if k == 0:
            ax.set_ylim(-1.35, 1.05)
            ax.legend(loc="upper right", fontsize=6.2, borderaxespad=0)
        else:
            ax.set_ylim(-0.35, 0.75)
    ax = fig.add_subplot(sub[0, 2])
    letter(ax, "c", x=-0.3, y=1.06)
    rows = []
    for comp, clab, col, mk in (("lncRNA TE-derived splice sites vs unused copies (main)", "unused copies", C["te"], "o"),
                                ("lncRNA TE-derived splice sites vs intronic", "intronic copies", C["te"], "s")):
        for cmp_, a_, lab, filled in ((comp, "ancestral=preformed", "preformed", True),
                                      (comp, "ancestral=created", "created", True),
                                      (comp + " (strata x CpG)", "ancestral=created", "created, CpG-matched", False)):
            q = S[(S.comparison == cmp_) & (S.subset == a_)]
            if len(q):
                r = q.iloc[0]
                rows.append((lab, r.mh_log_or, r.mh_log_or_lo, r.mh_log_or_hi, col, mk, filled))
    # rows: 3 subsets x 2 control sets -> grouped vertical bars
    labs = ["preformed", "created", "created,\nCpG strata"]
    for k, (lab, e, lo, hi, col, mk, filled) in enumerate(rows):
        gi, ci = k % 3, k // 3
        xx = gi + (ci - 0.5) * 0.36
        ax.bar(xx, e, width=0.34, color=C["te"] if ci == 0 else "#F0B48F", edgecolor=C["te"], lw=0.6)
        ax.plot([xx, xx], [lo, hi], color=INK, lw=0.8)
    ax.axhline(0, color=INK, lw=0.6)
    ax.set_xticks(range(3))
    ax.set_xticklabels(labs, fontsize=6.4)
    ax.set_xlim(-0.6, 2.6)
    ax.set_ylim(-0.4, 0.3)
    ax.set_ylabel("log OR of a variable site")
    ax.legend(handles=[Patch(facecolor=C["te"], edgecolor=C["te"], label="vs unused copies"),
                       Patch(facecolor="#F0B48F", edgecolor=C["te"], label="vs intronic copies")],
              loc="upper right", ncol=1, fontsize=6.0, borderaxespad=0.1, handlelength=1.2)
    ax.set_title("TE sites by origin")
    save(fig, "Fig7")


if __name__ == "__main__":
    from matplotlib.backends.backend_pdf import PdfPages
    want = set(sys.argv[4:]) if len(sys.argv) > 4 else None
    if want is None:
        PDF_ALL = PdfPages(os.path.join(OUT, "main_figures_v2.pdf"))
    for k, f in (("1", fig1), ("2", fig2), ("3", fig3), ("4", fig4), ("5", fig5), ("6", fig6)):
        if want is None or k in want:
            f()
    if PDF_ALL is not None:
        PDF_ALL.close()
