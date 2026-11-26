#!/usr/bin/env python3
"""Supplementary Figures S1-S7 and S9 of the manuscript (MBE numbering; S8 is drawn by make_figS8.py), Okabe-Ito colours.

Usage: python make_supp_figures.py <data_dir with work/> <revision table dir> <synthetic run dir>
                                      <shifted synthetic run dir> <out_dir> [figure numbers]
  synthetic run dir: tests/run_synthetic_test.sh output (work/prof_main, work/prof_pseudo)
  shifted run dir:   the same steps on synthetic data with a uniform 25% variant deficit in exonic sense Alu copies
                     (figures/make_synthetic_data_shift.py with SHIFT=0.75), work/prof_main
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np
import pandas as pd

D, RV, SYN, SHIFT, OUT = sys.argv[1:6]
WANT = set(sys.argv[6:])
W = os.path.join(D, "work")
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
OI = {"orange": "#E69F00", "sky": "#56B4E9", "green": "#009E73", "yellow": "#F0E442", "blue": "#0072B2",
      "vermillion": "#D55E00", "purple": "#CC79A7", "black": "#000000", "grey": "#999999"}
C = {"AluS": OI["orange"], "AluJ": OI["orange"], "AluY": OI["orange"], "MIR": OI["green"], "MIR3": OI["green"],
     "L2abc": OI["purple"], "exon": OI["vermillion"], "intron": "#8C8C8C", "acc": OI["vermillion"], "don": OI["blue"]}
LABEL = {"AluS": "AluS", "AluJ": "AluJ", "AluY": "AluY", "MIR": "MIR", "MIR3": "MIR3", "L2abc": "L2a/b/c"}
STAR = {"AluS": "AluS*", "AluJ": "AluJ*", "AluY": "AluY*", "MIR": "MIR*", "MIR3": "MIR3", "L2abc": "L2abc"}
MM = 1 / 25.4
PDF = PdfPages(os.path.join(OUT, "supplemental_figures_v2.pdf")) if not WANT else None


def letter(ax, s, x=-0.14, y=1.05):
    ax.text(x, y, s.upper(), transform=ax.transAxes, fontsize=10, fontweight="bold", va="bottom", color=INK)


def save(fig, name):
    fig.savefig(os.path.join(OUT, f"{name}.png"), dpi=300, facecolor="white")
    fig.savefig(os.path.join(OUT, f"{name}.pdf"), facecolor="white")
    if PDF is not None:
        PDF.savefig(fig, facecolor="white")
    plt.close(fig)


def rv(name):
    return pd.read_csv(os.path.join(RV, name), sep="\t")


def profile(ax, pr, col, label=None, centered=False, mark=True, lw=1.1):
    pr = pr.sort_values("bin")
    x = (pr["cons_from"] + pr["cons_to"]) / 2
    e, lo, hi = ("effect_centered", "cen_lo", "cen_hi") if centered else ("effect", "ci_lo", "ci_hi")
    ax.fill_between(x, pr[lo], pr[hi], color=col, alpha=0.2, lw=0)
    ax.plot(x, pr[e], color=col, lw=lw, label=label)
    if mark:
        q = pr["q_centered"] if centered else pr["q"]
        s = pr[q < 0.1]
        ax.plot((s["cons_from"] + s["cons_to"]) / 2, s[e], "o", ms=2.8, color=INK, zorder=4)
    ax.axhline(0, color=GRID, lw=0.8, zorder=0)


# ----------------------------------------------------------------------------- S1 simulations
def s1():
    fig = plt.figure(figsize=(180 * MM, 56 * MM))
    gs = GridSpec(1, 4, figure=fig, wspace=0.36, left=0.09, right=0.99, top=0.82, bottom=0.2)
    need = [os.path.join(d, "work", p, "profiles.tsv") for d in (SYN, SHIFT) for p in ("prof_main", "prof_pseudo")
            if d != "none"]
    miss = [p for p in need if not os.path.exists(p)]
    if SYN == "none" or miss:
        print("Fig. S1 skipped: synthetic runs of stage s10 not found (" + ", ".join(miss or [SYN]) + ")")
        return
    main = pd.read_csv(os.path.join(SYN, "work/prof_main/profiles.tsv"), sep="\t")
    pseudo = pd.read_csv(os.path.join(SYN, "work/prof_pseudo/profiles.tsv"), sep="\t")
    shift = pd.read_csv(os.path.join(SHIFT, "work/prof_main/profiles.tsv"), sep="\t")
    sel = lambda d, ori: d[(d.layer == "population") & (d.subfamily == "AluS*") & (d.orientation == ori)]
    ps = sel(pseudo, "sense")
    if not len(ps):
        ps = pseudo[(pseudo.layer == "population") & (pseudo.subfamily == "AluS*")]
    panels = [("a", "Exonic sense copies", sel(main, "sense"), C["exon"]),
              ("b", "Exonic antisense copies", sel(main, "antisense"), C["intron"]),
              ("c", "Pseudo-exons in introns", ps, C["intron"])]
    axes = []
    for i, (l, t, d, col) in enumerate(panels):
        ax = fig.add_subplot(gs[0, i], sharey=axes[0] if axes else None)
        axes.append(ax)
        ax.axvspan(101, 130, color="#FCE6CF", lw=0, zorder=0)
        profile(ax, d, col)
        ax.set_title(t)
        letter(ax, l)
        ax.set_xlabel("Alu consensus position")
        if i == 0:
            ax.set_ylabel("log OR of a variable site\n(below 0 = fewer variants)")
    ax = fig.add_subplot(gs[0, 3], sharey=axes[0])
    ax.axvspan(101, 130, color="#FCE6CF", lw=0, zorder=0)
    d = sel(shift, "sense")
    profile(ax, d, C["intron"], label="raw profile", mark=False)
    profile(ax, d, C["exon"], label="centered profile", centered=True)
    ax.set_title("Added uniform deficit")
    ax.set_xlabel("Alu consensus position")
    ax.legend(loc="lower right", fontsize=6.2, handlelength=1.2, borderaxespad=0.1)
    letter(ax, "d")
    lo = min(np.nanmin(x["ci_lo"]) for x in (sel(main, "sense"), d)) - 0.05
    axes[0].set_ylim(min(lo, -1.45), 0.75)
    save(fig, "FigS1")


# ----------------------------------------------------------------------------- S2 transcription starts
def s2():
    fig = plt.figure(figsize=(180 * MM, 64 * MM))
    gs = GridSpec(1, 2, figure=fig, wspace=0.22, left=0.1, right=0.8, top=0.85, bottom=0.22)
    regions = [("first exon", OI["vermillion"], "-"), ("internal exon", OI["orange"], "-"),
               ("last exon", OI["purple"], "-"), ("single exon", OI["black"], "--"),
               ("first intron", OI["blue"], "-"), ("other intron", OI["sky"], "-")]
    bins = ["0-250", "250-500", "500-1k", "1k-2k", "2k-5k", ">5k"]
    axes = []
    for i, (f, title) in enumerate((("tss_check.tsv", "Analyzed lncRNAs (random 4,000)"),
                                    ("tss_check_strict.tsv", "lncRNAs present in GTEx (random 6,000)"))):
        t = pd.read_csv(os.path.join(W, f), sep="\t")
        ax = fig.add_subplot(gs[0, i], sharey=axes[0] if axes else None)
        axes.append(ax)
        for reg, col, ls in regions:
            s = t[(t.region == reg) & (t.dist_from_tss.isin(bins)) & (t.bases >= 20000)].set_index("dist_from_tss")
            xs = [bins.index(b) for b in s.index]
            ax.plot(xs, s["mean_p447"], ls, marker="o", ms=3, color=col, lw=1.1, label=reg)
        ax.axhline(0, color=GRID, lw=0.8, zorder=0)
        ax.axvspan(-0.4, 1.5, color="#EEEEEE", lw=0, zorder=0)
        ax.set_xticks(range(len(bins)))
        ax.set_xticklabels([b.replace("-", "–") for b in bins])
        ax.set_xlim(-0.4, len(bins) - 0.6)
        ax.set_xlabel("distance from the transcription start (nt)")
        ax.set_title(title)
        letter(ax, "ab"[i])
        if i == 0:
            ax.set_ylabel("mean phyloP (447 mammals),\nrepeat-free sequence")
    axes[1].legend(loc="center left", bbox_to_anchor=(1.02, 0.5))
    save(fig, "FigS2")


# ----------------------------------------------------------------------------- S3 planted constraint, all settings
def s3():
    groups = [("AluS", "AluS*"), ("AluJ", "AluJ*"), ("MIR", "MIR*"), ("L2abc", "L2abc")]
    P = pd.concat([rv(f"power_{g}_{o}_power.tsv") for g, _ in groups for o in ("sense", "antisense")])
    reds = (0.0, 0.1, 0.25, 0.5)
    fracs = (0.25, 0.5, 1.0)
    fig = plt.figure(figsize=(180 * MM, 225 * MM))
    gs = GridSpec(8, 3, figure=fig, hspace=0.55, wspace=0.12, left=0.19, right=0.9, top=0.955, bottom=0.05)
    im = None
    for r, ((g, star), ori) in enumerate([(gg, o) for gg in groups for o in ("sense", "antisense")]):
        for c, phi in enumerate(fracs):
            ax = fig.add_subplot(gs[r, c])
            t = P[(P.group == star) & (P.orientation == ori)]
            segs = sorted(t.segment.unique(), key=lambda s_: (s_ == "all", 0 if s_ == "all" else
                                                              int(s_.split("-")[1]) - int(s_.split("-")[0])))
            M = np.full((len(segs), len(reds)), np.nan)
            for si, sg in enumerate(segs):
                for ci, rr in enumerate(reds):
                    q = t[(t.segment == sg) & (t.reduction == rr) & ((t.fraction_units == phi) if rr > 0 else True)]
                    if len(q):
                        M[si, ci] = q.power_whole.iloc[0] if sg == "all" else q.power_window.iloc[0]
            im = ax.imshow(M, cmap="cividis", vmin=0, vmax=1, aspect="auto")
            for si in range(M.shape[0]):
                for ci in range(M.shape[1]):
                    v = M[si, ci]
                    if np.isfinite(v):
                        ax.text(ci, si, f"{v:.2f}", ha="center",
                                va="center", fontsize=5.8, color="black" if v > 0.55 else "white")
            ax.axvline(0.5, color="white", lw=1.2)
            nts = [f"{int(sg.split('-')[1]) - int(sg.split('-')[0]) + 1} nt" if sg != "all" else "whole" for sg in segs]
            ax.set_yticks(range(len(segs)))
            ax.set_yticklabels(nts if c == 0 else [""] * len(segs), fontsize=6)
            ax.set_xticks(range(len(reds)))
            ax.set_xticklabels(["0"] + [f"−{int(100 * x)}%" for x in reds[1:]] if r == 7 else [""] * len(reds),
                               fontsize=6)
            for sp in ax.spines.values():
                sp.set_visible(False)
            ax.tick_params(length=0)
            if r == 0:
                ax.set_title(f"domain in {int(100 * phi)}% of exonic copies", fontsize=7, fontweight="normal",
                             loc="center", pad=8)
            if c == 0:
                ax.text(-0.3, 0.5, f"{LABEL[g]}\n{ori}", transform=ax.transAxes, ha="right", va="center",
                        fontsize=7, fontweight="bold", color=C[g])
            if r == 7 and c == 1:
                ax.set_xlabel("variants removed in the domain (0 = nothing planted: unmodified data)", fontsize=6.8)
    cax = fig.add_axes([0.925, 0.35, 0.012, 0.3])
    cb = fig.colorbar(im, cax=cax)
    cb.ax.tick_params(labelsize=6, length=1.5)
    cb.outline.set_linewidth(0.4)
    cb.set_label("detection rate (20 replicates)", fontsize=6.5)
    fig.text(0.015, 0.99, "A", fontsize=10, fontweight="bold", va="top")
    save(fig, "FigS3")


# ----------------------------------------------------------------------------- S4 profiles, all groups, same genes
def s4():
    groups = ["AluS", "AluJ", "AluY", "MIR", "MIR3", "L2abc"]
    cols = [("population", "sense"), ("population", "antisense"), ("phylop_p447", "sense"), ("phylop_p447", "antisense")]
    fig, axes = plt.subplots(len(groups), 4, figsize=(180 * MM, 215 * MM), sharex="row")
    fig.subplots_adjust(left=0.16, right=0.99, top=0.95, bottom=0.08, hspace=0.62, wspace=0.45)
    for r, g in enumerate(groups):
        base = os.path.join(W, "prof_v5", f"{g}_b20_paired")
        pr = pd.read_csv(os.path.join(base, "profiles.tsv"), sep="\t")
        pk = os.path.join(base, "peaks_centered.tsv")
        try:
            peaks = pd.read_csv(pk, sep="\t")
        except (pd.errors.EmptyDataError, FileNotFoundError):
            peaks = pd.DataFrame(columns=["layer", "orientation", "cons_from", "cons_to"])
        for c, (layer, ori) in enumerate(cols):
            ax = axes[r][c]
            d = pr[(pr.layer == layer) & (pr.orientation == ori)]
            for _, p in peaks[(peaks.layer == layer) & (peaks.orientation == ori)].iterrows():
                ax.axvspan(p.cons_from, p.cons_to, color="#DDDDDD", lw=0, zorder=0)
            if len(d) and d["effect"].notna().any():
                profile(ax, d.dropna(subset=["effect"]), C[g])
            else:
                ax.text(0.5, 0.5, "too few copies", transform=ax.transAxes, ha="center", fontsize=6, color="#777777")
            if r == 0:
                ax.set_title(("Variants, " if layer == "population" else "phyloP, ") + ori, fontsize=7.5)
            ax.set_ylabel(("log OR" if layer == "population" else "phyloP difference"), fontsize=6.5)
            if c == 0:
                bb = ax.get_position()
                fig.text(0.012, (bb.y0 + bb.y1) / 2, LABEL[g], ha="left", va="center", fontsize=7.5,
                         fontweight="bold", color=C[g])
            if r == len(groups) - 1:
                ax.set_xlabel("consensus position")
    fig.legend(handles=[Line2D([], [], color="#777777", lw=1.1, label="exon − intron, same genes (95% band)"),
                        Line2D([], [], marker="o", ls="", color=INK, ms=2.8, label="window differing from zero (q < 0.1)"),
                        Patch(color="#DDDDDD", label="segment differing from the rest of the element (q < 0.1)")],
               loc="lower center", ncol=3, bbox_to_anchor=(0.53, 0.0), fontsize=6.2)
    save(fig, "FigS4")


# ----------------------------------------------------------------------------- S5 subsets: whole-element estimates
def s5():
    groups = ["AluS", "AluJ", "MIR", "L2abc"]
    subsets = [("all", "all copies", "o", True), ("embedded", "embedded in an exon", "s", False),
               ("far50", "bases ≥ 50 nt from splice sites", "^", False), ("dimer", "complete Alu dimers", "D", False)]
    fig = plt.figure(figsize=(180 * MM, 95 * MM))
    gs = GridSpec(1, 2, figure=fig, wspace=0.12, left=0.13, right=0.98, top=0.86, bottom=0.2)
    rows = [(g, o) for g in groups for o in ("sense", "antisense")]
    y = np.arange(len(rows))[::-1]
    for k, (layer, xlab, title, xlim) in enumerate((
            ("population", "log OR of a variable site (below 0 = fewer variants)", "Variants", (-0.2, 0.15)),
            ("phylop_consmatch", "difference in the share of bases matching the consensus", "Consensus match",
             (-0.07, 0.05)))):
        ax = fig.add_subplot(gs[0, k])
        letter(ax, "ab"[k], x=-0.05 if k else -0.28)
        ax.axvline(0, color=GRID, lw=0.8, zorder=0)
        for yi, (g, o) in zip(y, rows):
            for j, (sub, lab, mk, filled) in enumerate(subsets):
                f = os.path.join(RV, f"whole_prof_sub_{g}_{sub}.tsv")
                if not os.path.exists(f):
                    continue
                w = pd.read_csv(f, sep="\t")
                w = w[(w.layer == layer) & (w.orientation == o)]
                if not len(w):
                    continue
                r = w.iloc[0]
                yy = yi + 0.3 - j * 0.2
                ax.plot([r.ci_lo, r.ci_hi], [yy] * 2, color=C[g], lw=1.0)
                ax.plot(r.effect, yy, mk, ms=3.6, color=C[g], mfc=C[g] if filled else "white", mew=0.9)
        ax.set_yticks(y)
        ax.set_yticklabels([f"{LABEL[g]} {o}" for g, o in rows] if k == 0 else [""] * len(rows))
        ax.tick_params(axis="y", length=0)
        ax.set_ylim(-0.7, len(rows) - 0.3)
        ax.set_xlim(*xlim)
        ax.set_xlabel(xlab)
        ax.set_title(title)
    fig.legend(handles=[Line2D([], [], marker=mk, ls="", color="#555555", mfc="#555555" if f else "white", ms=3.8,
                               label=lab) for _, lab, mk, f in subsets],
               loc="upper center", ncol=4, bbox_to_anchor=(0.55, 1.0), fontsize=6.5)
    save(fig, "FigS5")


# ----------------------------------------------------------------------------- S6 conserved copies by compartment
def s6():
    t = rv("conserved_p447_paired.tsv")
    fams = ["MIR", "MIRb", "MIRc", "MIR3", "MIR1_Amn", "L2", "L2a", "L2b", "L2c", "L2d", "L2d2"]
    fig = plt.figure(figsize=(180 * MM, 80 * MM))
    gs = GridSpec(1, 2, figure=fig, wspace=0.12, left=0.12, right=0.98, top=0.85, bottom=0.18)
    y = np.arange(len(fams))[::-1]
    for k, des in enumerate(("all genes", "same genes")):
        ax = fig.add_subplot(gs[0, k])
        letter(ax, "ab"[k], x=-0.05 if k else -0.2)
        for yi, f in zip(y, fams):
            r = t[(t.family == f) & (t.design == des)]
            if not len(r):
                continue
            r = r.iloc[0]
            col = C["MIR"] if f.startswith("MIR") else C["L2abc"]
            ax.plot([100 * r.share_intron, 100 * r.share_exon], [yi, yi], color=GRID, lw=1.2, zorder=1)
            ax.plot(100 * r.share_intron, yi, "s", ms=3.6, color=col, mfc="white", mew=0.9, zorder=2)
            ax.plot(100 * r.share_exon, yi, "o", ms=3.8, color=col, zorder=3)
            ortxt = (f"OR {r.mh_or:.2f} [{r.ci_lo:.2f}–{r.ci_hi:.2f}]" if np.isfinite(r.ci_hi) else "not testable")
            ax.text(10.4, yi, f"{int(r.conserved_exon)}/{int(r.units_exon)} · {int(r.conserved_intron)}/{int(r.units_intron)}   {ortxt}",
                    va="center", fontsize=5.6, color="#555555")
        ax.set_yticks(y)
        ax.set_yticklabels(fams if k == 0 else [""] * len(fams))
        ax.tick_params(axis="y", length=0)
        ax.set_xlim(0, 21)
        ax.spines['bottom'].set_bounds(0, 10)
        ax.set_xticks([0, 2, 4, 6, 8, 10])
        ax.set_ylim(-0.7, len(fams) - 0.3)
        ax.set_xlabel("units with a conserved 20-nt stretch (%)")
        ax.set_title("All genes" if k == 0 else "Same genes")
    fig.legend(handles=[Line2D([], [], marker="o", ls="", color="#555555", ms=3.8, label="exonic copies"),
                        Line2D([], [], marker="s", ls="", mfc="white", mec="#555555", ms=3.6, label="intronic copies"),
                        Line2D([], [], ls="", label="numbers: conserved / units (exonic · intronic); odds ratio, GC- and length-stratified")],
               loc="upper center", ncol=3, bbox_to_anchor=(0.55, 1.0), fontsize=6.5)
    save(fig, "FigS6")


# ----------------------------------------------------------------------------- S7 RNA properties of host genes
def s7():
    T = pd.read_csv(os.path.join(W, "trend", "te_vs_nonte_trend.tsv"), sep="\t")
    props = [("stab", "RNA stability"), ("nuc", "nuclear retention"), ("splice", "splicing efficiency"),
             ("synth", "synthesis rate"), ("gtex", "expression (GTEx)"), ("tau", "tissue specificity (tau)")]
    groups = ["AluS", "AluJ", "MIR", "L2abc"]
    fig = plt.figure(figsize=(180 * MM, 78 * MM))
    gs = GridSpec(1, 2, figure=fig, wspace=0.1, left=0.15, right=0.98, top=0.85, bottom=0.17)
    y = np.arange(len(props))[::-1]
    for k, (layer, title) in enumerate((("excess_poly_noncpg", "TE-specific excess, variants"),
                                        ("excess_phylop_p447", "TE-specific excess, phyloP (mammals)"))):
        ax = fig.add_subplot(gs[0, k])
        letter(ax, "ab"[k], x=-0.05 if k else -0.35)
        ax.axvline(0, color=GRID, lw=0.8, zorder=0)
        for j, g in enumerate(groups):
            for yi, (p, _) in zip(y, props):
                r = T[(T.subfamily == STAR[g]) & (T.layer == layer) & (T.orientation == "pooled") & (T.property == p)]
                if not len(r):
                    continue
                r = r.iloc[0]
                yy = yi + 0.27 - j * 0.18
                ax.plot([r.ci_lo, r.ci_hi], [yy] * 2, color=C[g], lw=1.0)
                ax.plot(r.spearman_rho, yy, ["o", "s", "o", "o"][j], ms=3.4, color=C[g],
                        mfc=C[g] if j != 1 else "white", mew=0.9)
        ax.set_yticks(y)
        ax.set_yticklabels([lab for _, lab in props] if k == 0 else [""] * len(props))
        ax.tick_params(axis="y", length=0)
        ax.set_xlim(-0.65, 0.65)
        ax.set_ylim(-0.7, len(props) - 0.3)
        ax.set_xlabel("Spearman correlation with the host-gene property")
        ax.set_title(title)
    fig.legend(handles=[Line2D([], [], marker="o" if g != "AluJ" else "s", ls="", color=C[g],
                               mfc=C[g] if g != "AluJ" else "white", ms=3.6, label=LABEL[g]) for g in groups],
               loc="upper center", ncol=4, bbox_to_anchor=(0.56, 1.0), fontsize=6.5)
    save(fig, "FigS7")


# ----------------------------------------------------------------------------- S8 MIR and L2 splice positions
def s8():
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
            for feat, sign, col in (("acceptor", 1, C["acc"]), ("donor", -1, C["don"])):
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
            ax.set_ylabel("lncRNA splice sites")
            ax.set_title(f"{fam} copies {ori} to the transcript")
            if r == 0 and c == 1:
                ax.legend(handles=[Patch(color=C["acc"], label="acceptors (up)"),
                                   Patch(color=C["don"], label="donors (down)")], loc="lower right", handlelength=1)
    save(fig, "FigS9")  # cited after the LTR figure (Fig. S8, make_figS8.py)


if __name__ == "__main__":
    for k, f in (("1", s1), ("2", s2), ("3", s3), ("4", s4), ("5", s5), ("6", s6), ("7", s7), ("8", s8)):
        if not WANT or k in WANT:
            f()
    if PDF is not None:
        PDF.close()
