#!/usr/bin/env python3
"""Supplementary figures S1-S4 for the RIDL 2.0 manuscript.

Usage: python make_supp_figures.py <data_dir with work/...> <synthetic run dir> <shifted synthetic run dir> <out_dir>
  synthetic run dir: output of tests/run_synthetic_test.sh (work/prof_main, work/prof_pseudo)
  shifted run dir:   the same pipeline on synthetic data with a uniform 25% variant deficit in exonic sense Alu
                     copies (tests/make_synthetic_data.py with SHIFT=0.75), step 6 output in work/prof_main
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np
import pandas as pd

D, SYN, SHIFT, OUT = sys.argv[1:5]
W = os.path.join(D, "work")
os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 7, "axes.titlesize": 7.5, "axes.labelsize": 7,
    "xtick.labelsize": 6.5, "ytick.labelsize": 6.5, "legend.fontsize": 6.5, "axes.linewidth": 0.6,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6, "xtick.major.size": 2.5, "ytick.major.size": 2.5,
    "axes.spines.top": False, "axes.spines.right": False, "pdf.fonttype": 42,
    "axes.titleweight": "bold", "axes.titlelocation": "left", "legend.frameon": False,
})
INK, GRID = "#222222", "#d9d9d9"
C = {"AluS": "#E69F00", "AluJ": "#F2C16B", "AluY": "#A86F00", "MIR": "#009E73", "MIR3": "#66C2A5",
     "L2abc": "#CC79A7", "exon": "#A50F15", "intron": "#8C8C8C"}
LABEL = {"AluS": "Alu S", "AluJ": "Alu J", "AluY": "Alu Y", "MIR": "MIR", "MIR3": "MIR3", "L2abc": "L2 a/b/c"}
MM = 1 / 25.4
PDF = PdfPages(os.path.join(OUT, "RIDL2_supplementary_figures.pdf"))


def letter(ax, s, x=-0.14, y=1.05):
    ax.text(x, y, s, transform=ax.transAxes, fontsize=10, fontweight="bold", va="bottom", color=INK)


def save(fig, name):
    fig.savefig(os.path.join(OUT, f"{name}.png"), dpi=300, facecolor="white")
    fig.savefig(os.path.join(OUT, f"{name}.pdf"), facecolor="white")
    PDF.savefig(fig, facecolor="white")
    plt.close(fig)


def profile(ax, pr, col, label=None, centered=False, mark=True):
    pr = pr.sort_values("bin")
    x = (pr["cons_from"] + pr["cons_to"]) / 2
    e, lo, hi = ("effect_centered", "cen_lo", "cen_hi") if centered else ("effect", "ci_lo", "ci_hi")
    ax.fill_between(x, pr[lo], pr[hi], color=col, alpha=0.18, lw=0)
    ax.plot(x, pr[e], color=col, lw=1.1, label=label)
    if mark:
        q = pr["q_centered"] if centered else pr["q"]
        s = pr[q < 0.1]
        ax.plot((s["cons_from"] + s["cons_to"]) / 2, s[e], "o", ms=2.8, color=INK, zorder=4)
    ax.axhline(0, color=GRID, lw=0.8, zorder=0)


# ----------------------------------------------------------------------------- S1
def s1():
    fig = plt.figure(figsize=(180 * MM, 52 * MM))
    gs = GridSpec(1, 4, figure=fig, wspace=0.36, left=0.09, right=0.99, top=0.82, bottom=0.2)
    main = pd.read_csv(os.path.join(SYN, "work/prof_main/profiles.tsv"), sep="\t")
    pseudo = pd.read_csv(os.path.join(SYN, "work/prof_pseudo/profiles.tsv"), sep="\t")
    shift = pd.read_csv(os.path.join(SHIFT, "work/prof_main/profiles.tsv"), sep="\t")
    sel = lambda d, ori: d[(d.layer == "population") & (d.subfamily == "AluS*") & (d.orientation == ori)]
    panels = [("a", "Sense copies", sel(main, "sense")),
              ("b", "Antisense copies", sel(main, "antisense")),
              ("c", "Pseudo-exons", sel(pseudo, "sense") if len(sel(pseudo, "sense")) else
               pseudo[(pseudo.layer == "population") & (pseudo.subfamily == "AluS*")])]
    axes = []
    for i, (l, t, d) in enumerate(panels):
        ax = fig.add_subplot(gs[0, i], sharey=axes[0] if axes else None)
        axes.append(ax)
        ax.axvspan(101, 130, color="#FCE6CF", lw=0, zorder=0)
        profile(ax, d, C["exon"] if i == 0 else C["intron"])
        ax.set_title(t)
        letter(ax, l)
        ax.set_xlabel("Alu consensus position")
        if i == 0:
            ax.set_ylabel("log OR of a variable site\n(below 0 = fewer variants in exons)")
    ax = fig.add_subplot(gs[0, 3], sharey=axes[0])
    ax.axvspan(101, 130, color="#FCE6CF", lw=0, zorder=0)
    d = sel(shift, "sense")
    profile(ax, d, C["intron"], label="raw (shifted)", mark=False)
    profile(ax, d, C["exon"], label="centered", centered=True)
    ax.set_title("Shift + segment")
    ax.set_xlabel("Alu consensus position")
    ax.legend(loc="lower left", fontsize=5.6, handlelength=1.2, borderaxespad=0.1)
    letter(ax, "d")
    axes[0].set_ylim(-1.45, 0.45)
    save(fig, "FigS1_simulation")


# ----------------------------------------------------------------------------- S2
def s2():
    fig = plt.figure(figsize=(180 * MM, 62 * MM))
    gs = GridSpec(1, 2, figure=fig, wspace=0.22, left=0.1, right=0.8, top=0.85, bottom=0.22)
    regions = [("first exon", "#A50F15", "-"), ("internal exon", "#E3646A", "-"), ("last exon", "#F2A3A6", "-"),
               ("single exon", "#7A0A10", "--"), ("first intron", "#555555", "-"), ("other intron", "#A6A6A6", "-")]
    bins = ["0-250", "250-500", "500-1k", "1k-2k", "2k-5k", ">5k"]
    axes = []
    for i, (f, title) in enumerate((("tss_check.tsv", "All analysed lncRNAs (4,000 genes)"),
                                    ("tss_check_strict.tsv", "Well-expressed lncRNAs (6,000 genes)"))):
        t = pd.read_csv(os.path.join(W, f), sep="\t")
        ax = fig.add_subplot(gs[0, i], sharey=axes[0] if axes else None)
        axes.append(ax)
        for reg, col, ls in regions:
            s = t[(t.region == reg) & (t.dist_from_tss.isin(bins)) & (t.bases >= 20000)].set_index("dist_from_tss")
            xs = [bins.index(b) for b in s.index]
            ax.plot(xs, s["mean_p447"], ls, marker="o", ms=3, color=col, lw=1.1, label=reg)
        ax.axhline(0, color=GRID, lw=0.8, zorder=0)
        ax.set_xticks(range(len(bins)))
        ax.set_xticklabels([b.replace("-", "–") for b in bins])
        ax.set_xlabel("distance from the transcription start (nt)")
        ax.set_title(title)
        letter(ax, "ab"[i])
        if i == 0:
            ax.set_ylabel("mean phyloP (447 mammals),\nrepeat-free sequence")
    axes[1].legend(loc="center left", bbox_to_anchor=(1.02, 0.5))
    save(fig, "FigS2_tss")


# ----------------------------------------------------------------------------- S3
def s3():
    fig = plt.figure(figsize=(180 * MM, 105 * MM))
    gs = GridSpec(2, 2, figure=fig, wspace=0.2, hspace=0.55, left=0.08, right=0.99, top=0.93, bottom=0.09)
    hs = pd.read_csv(os.path.join(W, "side/architecture_hotspots.tsv"), sep="\t")
    hs = hs[hs.gene_set == "lnc"]
    for r, (fam, binw, xlim) in enumerate((("MIR", 3, (0, 270)), ("L2", 20, (2000, 3460)))):
        for c, ori in enumerate(("antisense", "sense")):
            ax = fig.add_subplot(gs[r, c])
            if c == 0:
                letter(ax, "ab"[r], x=-0.12)
            bins = np.arange(xlim[0], xlim[1] + binw, binw)
            for feat, sign, col in (("acceptor", 1, "#C0392B"), ("donor", -1, "#2E6DA4")):
                d = hs[(hs.rfamily == fam) & (hs.orientation == ori) & (hs.feature == feat)]
                h, _ = np.histogram(d.cons_pos, bins=bins, weights=d.n)
                ax.bar(bins[:-1] + binw / 2, sign * h, width=binw, color=col, lw=0)
            ax.axhline(0, color=INK, lw=0.5)
            ax.set_xlim(*xlim)
            yl = ax.get_ylim()
            ax.set_yticks(ax.get_yticks())
            ax.set_yticklabels([f"{abs(int(t))}" for t in ax.get_yticks()])
            ax.set_ylim(yl)
            ax.set_xlabel(f"{fam} consensus position" + (" (3′ part)" if fam == "L2" else ""))
            if c == 0:
                ax.set_ylabel("lncRNA splice sites")
            ax.set_title(f"{fam} copies {ori} to the transcript")
            if r == 0 and c == 1:
                ax.legend(handles=[Patch(color="#C0392B", label="acceptors"), Patch(color="#2E6DA4", label="donors")],
                          loc="lower right", handlelength=1)
    save(fig, "FigS3_mir_l2_splice")


# ----------------------------------------------------------------------------- S4
def s4():
    groups = ["AluS", "AluJ", "AluY", "MIR", "MIR3", "L2abc"]
    cols = [("population", "sense"), ("population", "antisense"), ("phylop_p447", "sense"), ("phylop_p447", "antisense")]
    fig, axes = plt.subplots(len(groups), 4, figsize=(180 * MM, 215 * MM), sharex="row")
    fig.subplots_adjust(left=0.1, right=0.99, top=0.95, bottom=0.075, hspace=0.6, wspace=0.42)
    for r, g in enumerate(groups):
        pr = pd.read_csv(os.path.join(W, "prof_v5", f"{g}_b20_paired", "profiles.tsv"), sep="\t")
        for c, (layer, ori) in enumerate(cols):
            ax = axes[r][c]
            d = pr[(pr.layer == layer) & (pr.orientation == ori)]
            if len(d) and d["effect"].notna().any():
                profile(ax, d.dropna(subset=["effect"]), C[g])
                s = d[d["q_centered"] < 0.1]
                ax.plot((s["cons_from"] + s["cons_to"]) / 2, s["effect"], "o", ms=3.2, mfc="none", mec=INK, mew=0.8)
            else:
                ax.text(0.5, 0.5, "too few copies", transform=ax.transAxes, ha="center", fontsize=6, color="#777777")
            if r == 0:
                ax.set_title(("Variants, " if layer == "population" else "phyloP, ") + ori, fontsize=7)
            if c == 0:
                ax.set_ylabel(f"{LABEL[g]}\n" + "log OR", fontsize=6.8)
            if c == 2:
                ax.set_ylabel("phyloP difference", fontsize=6.3)
            if r == len(groups) - 1:
                ax.set_xlabel("consensus position")
    fig.legend(handles=[Line2D([], [], color="#777777", lw=1.1, label="exon − intron, same genes (95% band)"),
                        Line2D([], [], marker="o", ls="", color=INK, ms=2.8, label="differing from zero (q < 0.1)"),
                        Line2D([], [], marker="o", ls="", mfc="none", mec=INK, ms=3.2,
                               label="differing from the rest of the element (q < 0.1)")],
               loc="lower center", ncol=3, bbox_to_anchor=(0.53, -0.005), fontsize=6.0)
    save(fig, "FigS4_profiles_all_groups")


if __name__ == "__main__":
    s1()
    s2()
    s3()
    s4()
    PDF.close()
