#!/usr/bin/env python3
"""Main figures for the RIDL 2.0 manuscript, drawn from the result tables of the pipeline.

Usage: python make_figures.py <data_dir with work/...> <out_dir> [fig numbers...]
Each figure is written as PNG (300 dpi) and PDF (vector).
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrow, Patch, Rectangle
import numpy as np
import pandas as pd

D = sys.argv[1] if len(sys.argv) > 1 else "data"
OUT = sys.argv[2] if len(sys.argv) > 2 else "out"
W = os.path.join(D, "work")
os.makedirs(OUT, exist_ok=True)

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 7, "axes.titlesize": 7.5, "axes.labelsize": 7,
    "xtick.labelsize": 6.5, "ytick.labelsize": 6.5, "legend.fontsize": 6.5, "axes.linewidth": 0.6,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6, "xtick.major.size": 2.5, "ytick.major.size": 2.5,
    "axes.spines.top": False, "axes.spines.right": False, "pdf.fonttype": 42, "svg.fonttype": "none",
    "axes.titleweight": "bold", "axes.titlelocation": "left", "legend.frameon": False,
})
INK = "#222222"
GRID = "#d9d9d9"
C = {"AluS": "#E69F00", "AluJ": "#F2C16B", "AluY": "#A86F00", "MIR": "#009E73", "MIR3": "#66C2A5",
     "L2abc": "#CC79A7", "ord": "#4D4D4D", "pc": "#0072B2", "exon": "#A50F15", "intron": "#8C8C8C",
     "intergenic": "#C9C9C9"}
LABEL = {"AluS": "Alu S", "AluJ": "Alu J", "AluY": "Alu Y", "MIR": "MIR", "MIR3": "MIR3", "L2abc": "L2 a/b/c"}
GROUPS = ["AluS", "AluJ", "AluY", "MIR", "MIR3", "L2abc"]
MM = 1 / 25.4


def tsv(*p):
    return pd.read_csv(os.path.join(W, *p), sep="\t")


def letter(ax, s, x=-0.12, y=1.04):
    ax.text(x, y, s, transform=ax.transAxes, fontsize=10, fontweight="bold", va="bottom", ha="left", color=INK)


PDF_ALL = None


def save(fig, name):
    fig.savefig(os.path.join(OUT, f"{name}.png"), dpi=300, facecolor="white")
    fig.savefig(os.path.join(OUT, f"{name}.pdf"), facecolor="white")
    if PDF_ALL is not None:
        PDF_ALL.savefig(fig, facecolor="white")
    plt.close(fig)


def whole(run, layer="population"):
    w = tsv("prof_v5", run, "whole_element.tsv")
    return w[w.layer == layer].set_index("orientation")


def ladder(ax, rows, xlim, xlabel, ref=0.0):
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
            ax.text(start - 0.02 * span if right else start + 0.02 * span, yi,
                    f"{e:+.2f} [{lo:.2f}, {hi:.2f}]", ha="right" if right else "left", va="center",
                    fontsize=6, color=col)
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
    fig = plt.figure(figsize=(180 * MM, 92 * MM))
    gs = GridSpec(2, 2, figure=fig, width_ratios=[1.35, 1], height_ratios=[1, 1], wspace=0.42, hspace=0.75,
                  left=0.03, right=0.98, top=0.93, bottom=0.11)
    # a: schematic
    ax = fig.add_subplot(gs[:, 0])
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    ax.axis("off")
    letter(ax, "a", x=0.0, y=1.0)
    # gene model
    ax.plot([4, 96], [90, 90], color=C["intron"], lw=1.2, zorder=1)
    for x0, x1 in ((4, 20), (42, 60), (82, 96)):
        ax.add_patch(Rectangle((x0, 86.5), x1 - x0, 7, color=C["exon"], zorder=2))

    def te(x, y, w, right=True, col=C["AluS"]):
        if right:
            ax.add_patch(FancyArrow(x, y, w, 0, width=2.6, head_width=4.2, head_length=2.2,
                                    length_includes_head=True, color=col, zorder=3, lw=0))
        else:
            ax.add_patch(FancyArrow(x + w, y, -w, 0, width=2.6, head_width=4.2, head_length=2.2,
                                    length_includes_head=True, color=col, zorder=3, lw=0))
    te(7, 90, 9, right=False)
    te(45, 90, 9)
    te(85, 90, 8, right=False)
    te(25, 90, 7)
    te(33, 90, 6, right=False)
    te(66, 90, 8)
    ax.text(89, 81, "exon: mature RNA", ha="center", fontsize=6.3, color=C["exon"])
    ax.text(30, 81, "intron: spliced out", ha="center", fontsize=6.3, color="#666666")
    ax.text(4, 99, "lncRNA gene with TE copies", fontsize=6.8, fontweight="bold", color=INK)
    # copies on consensus coordinates
    ax.text(4, 73, "Copies placed on consensus coordinates", fontsize=6.8, fontweight="bold", color=INK)
    rng = np.random.default_rng(4)
    x0, x1 = 14, 92
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
    ax.text(12, 62.5, "exonic", ha="right", va="center", fontsize=6.3, color=C["exon"])
    ax.text(12, 47, "intronic", ha="right", va="center", fontsize=6.3, color="#666666")
    ax.plot([x0, x1], [38.5, 38.5], color=INK, lw=0.6)
    for p in (1, 100, 200, 300):
        ax.plot([cx(p)] * 2, [38.5, 37], color=INK, lw=0.6)
        ax.text(cx(p), 34, str(p), ha="center", fontsize=5.8)
    ax.text(x1, 30.5, "consensus position", ha="right", fontsize=5.8, color="#555555")
    ax.plot(cx(252), 73.6, "o", ms=1.6, color=INK)
    ax.text(cx(252) + 1.5, 73.6, "variant", fontsize=5.8, va="center", color="#555555")
    # profile
    ax.text(4, 25, "Compare, position by position", fontsize=6.8, fontweight="bold", color=INK)
    p = np.arange(1, 301, 5)
    base = 13 + 0.6 * np.sin(p * 0.37)
    exo = base.copy()
    dip = np.exp(-0.5 * ((p - 140) / 9) ** 2)
    exo = 13 + 0.6 * np.sin(p * 0.29 + 1) - 7 * dip
    ax.add_patch(Rectangle((cx(win[0]), 3), cx(win[1]) - cx(win[0]), 18, color="#FCE6CF", zorder=0, lw=0))
    ax.plot(cx(p), base, color=C["intron"], lw=1.2)
    ax.plot(cx(p), exo, color=C["exon"], lw=1.4)
    ax.text(cx(win[1]) + 1.5, 6, "fewer variants in exonic copies\n= candidate RNA domain", fontsize=5.8,
            va="center", color=C["exon"])
    ax.text(x0 - 2, 13, "variant\nfrequency", ha="right", va="center", fontsize=5.8, color="#555555")

    # b: exonic units per group (all genes)
    axb = fig.add_subplot(gs[0, 1])
    letter(axb, "b", x=-0.30)
    rows = []
    for g in GROUPS:
        w = whole(f"{g}_b20")
        rows.append((g, int(w.loc["sense", "units_fg"]), int(w.loc["antisense", "units_fg"])))
    y = np.arange(len(rows))[::-1]
    for yi, (g, s, a) in zip(y, rows):
        axb.barh(yi + 0.19, s, height=0.36, color=C[g])
        axb.barh(yi - 0.19, a, height=0.36, color=C[g], alpha=0.45)
        axb.text(max(s, a) + 20, yi, f"{s:,} / {a:,}", va="center", fontsize=5.8)
    axb.set_yticks(y)
    axb.set_yticklabels([LABEL[r[0]] for r in rows])
    axb.tick_params(axis="y", length=0)
    axb.set_xlabel("exonic TE units (copies or their exonic parts)")
    axb.set_xlim(0, 1350)
    axb.legend(handles=[Patch(color="#999999", label="sense"),
                        Patch(color="#999999", alpha=0.45, label="antisense")],
               loc="lower right", ncol=1, handlelength=1, borderaxespad=0)
    axb.set_title("Material: 19,040 lncRNA loci")

    # c: detectable effect vs window size
    axc = fig.add_subplot(gs[1, 1])
    letter(axc, "c", x=-0.30)
    sizes = {10: "power", 20: "power_b20", 30: "power_b30"}
    key = {"AluS": "AluS*", "AluJ": "AluJ*", "AluY": "AluY*", "MIR": "MIR*", "MIR3": "MIR3", "L2abc": "L2abc"}
    for g in GROUPS:
        ys = []
        for s, d in sizes.items():
            t = tsv(d, "power_summary.tsv")
            t = t[t.subfamily == key[g]]
            ys.append(100 * t["median_rel_mde"].mean())
        axc.plot(list(sizes), ys, "-o", color=C[g], ms=3, lw=1.1, label=LABEL[g])
    axc.axhline(10, color=INK, lw=0.6, ls="--")
    axc.text(30.6, 10.3, "10%", fontsize=5.8, va="bottom", ha="left")
    axc.set_xticks([10, 20, 30])
    axc.set_xlim(8, 34)
    axc.set_ylim(0, 31)
    axc.set_xlabel("window along the consensus (nt)")
    axc.set_ylabel("smallest detectable reduction\nin exonic variants (%)")
    axc.legend(ncol=2, loc="upper right", handlelength=1.2, columnspacing=0.8, borderaxespad=0)
    axc.set_title("Power per window")
    save(fig, "Fig1_approach")


# ----------------------------------------------------------------------------- Figure 2
def fig2():
    fig = plt.figure(figsize=(180 * MM, 78 * MM))
    gs = GridSpec(1, 3, figure=fig, width_ratios=[0.8, 1.25, 1.0], wspace=0.62, left=0.085, right=0.985,
                  top=0.86, bottom=0.25)
    # a: gradient intergenic -> intron -> exon (variants, all genes)
    ax = fig.add_subplot(gs[0, 0])
    letter(ax, "a", x=-0.45)
    nig = tsv("prof_nonte", "intron_vs_intergenic", "whole_element.tsv")
    nex = tsv("prof_nonte", "exon_vs_intron", "whole_element.tsv")
    series = [("ordinary sequence", C["ord"],
               nig.query("layer=='population'").effect.mean(), nex.query("layer=='population'").effect.mean())]
    for g in ("AluS", "MIR"):
        ig = tsv("diag6", f"{g}_intron_vs_intergenic", "whole_element.tsv").query("layer=='population'").effect.mean()
        ex = whole(f"{g}_b20").effect.mean()
        series.append((LABEL[g], C[g], ig, ex))
    for lab, col, ig, ex in series:
        ys = [0, ig, ig + ex]
        ax.plot([0, 1, 2], ys, "-o", color=col, ms=3.5, lw=1.2, label=lab)
    ax.axhline(0, color=GRID, lw=0.8, zorder=0)
    ax.set_xticks([0, 1, 2])
    ax.set_xticklabels(["intergenic", "intron", "exon"])
    ax.set_xlim(-0.3, 2.3)
    ax.set_ylabel("variant odds relative to\nintergenic copies (log OR)")
    ax.legend(loc="lower left", handlelength=1.4, borderaxespad=0)
    ax.set_title("Variation falls with\ngenomic location")

    # b: across genes vs within genes
    axb = fig.add_subplot(gs[0, 1])
    letter(axb, "b", x=-0.42)
    rows = []
    for g in GROUPS:
        a, p = whole(f"{g}_b20"), whole(f"{g}_b20_paired")
        for ori, tag in (("sense", "s"), ("antisense", "as")):
            rows.append((f"{LABEL[g]} {tag}", a.loc[ori], p.loc[ori], C[g]))
    o_all = nex.query("layer=='population'").iloc[0]
    o_par = tsv("prof_nonte", "exon_vs_intron_paired", "whole_element.tsv").query("layer=='population'").iloc[0]
    rows.append(("ordinary seq.", o_all, o_par, C["ord"]))
    n = len(rows)
    y = np.arange(n)[::-1]
    axb.axvspan(o_par.ci_lo, o_par.ci_hi, color="#E8E8E8", lw=0, zorder=0)
    axb.axvline(0, color=GRID, lw=0.8, zorder=0)
    for yi, (lab, A, P, col) in zip(y, rows):
        axb.plot([A.ci_lo, A.ci_hi], [yi + 0.18] * 2, color=col, lw=1.0, alpha=0.55)
        axb.plot(A.effect, yi + 0.18, "o", ms=3.4, mfc="white", mec=col, mew=0.9)
        axb.plot([P.ci_lo, P.ci_hi], [yi - 0.18] * 2, color=col, lw=1.2)
        axb.plot(P.effect, yi - 0.18, "o", ms=3.6, color=col)
    axb.set_yticks(y)
    axb.set_yticklabels([r[0] for r in rows], fontsize=6.2)
    axb.tick_params(axis="y", length=0)
    axb.set_xlim(-0.22, 0.28)
    axb.set_ylim(-0.7, n - 0.3)
    axb.set_xlabel("exon vs intron, log OR of a variable site")
    axb.legend(handles=[Line2D([], [], marker="o", ls="", mfc="white", mec="#555555", ms=3.4, label="all genes"),
                        Line2D([], [], marker="o", ls="", color="#555555", ms=3.6, label="same genes"),
                        Patch(color="#E8E8E8", label="ordinary seq., same genes")],
               loc="upper center", bbox_to_anchor=(0.42, -0.21), ncol=3, fontsize=5.8, handletextpad=0.3,
               columnspacing=0.9, borderaxespad=0)
    axb.set_title("Within genes, the deficit vanishes")

    # c: phyloP by GC bin, exon vs intron (ordinary sequence, same genes)
    sub = GridSpecFromSubplotSpec(2, 1, subplot_spec=gs[0, 2], height_ratios=[2.2, 1], hspace=0.12)
    axc = fig.add_subplot(sub[0])
    axd = fig.add_subplot(sub[1], sharex=axc)
    letter(axc, "c", x=-0.32)
    b = tsv("comp", "lnc_paired_bins.tsv")
    x = np.arange(len(b))
    axc.plot(x, b.phylop_p447_exon, "-o", color=C["exon"], ms=3.4, lw=1.2, label="exon")
    axc.plot(x, b.phylop_p447_intron, "-o", color=C["intron"], ms=3.4, lw=1.2, label="intron")
    axc.set_ylabel("mean phyloP")
    axc.legend(loc="upper right", borderaxespad=0)
    s = tsv("comp", "lnc_paired_summary.tsv").set_index(["metric", "comparison"])
    um, gm = s.loc[("phylop_p447", "unmatched")], s.loc[("phylop_p447", "gc_matched")]
    axc.text(0.03, 0.04, f"exon − intron\nunmatched: {um.exon_minus_intron:+.3f}\n"
             f"GC-matched: {gm.exon_minus_intron:+.3f}\n[{gm.ci_lo:.3f}, {gm.ci_hi:.3f}]",
             transform=axc.transAxes, fontsize=5.8, va="bottom")
    plt.setp(axc.get_xticklabels(), visible=False)
    axc.set_title("GC content explains lower\nexonic conservation")
    fe = b.units_exon / b.units_exon.sum()
    fi = b.units_intron / b.units_intron.sum()
    axd.bar(x - 0.18, fe, width=0.36, color=C["exon"])
    axd.bar(x + 0.18, fi, width=0.36, color=C["intron"])
    axd.set_ylabel("share of\nwindows", fontsize=6.3)
    axd.set_xticks(x)
    axd.set_xticklabels([r.replace("-", "–") for r in b.gc_range], fontsize=5.6, rotation=35, ha="right")
    axd.set_xlabel("GC content of 100-nt windows")
    save(fig, "Fig2_confounders")



# ----------------------------------------------------------------------------- Figure 3
def fig3():
    fig = plt.figure(figsize=(180 * MM, 98 * MM))
    gs = GridSpec(2, 2, figure=fig, width_ratios=[1, 1], height_ratios=[1, 1.05], wspace=0.55, hspace=0.62,
                  left=0.12, right=0.92, top=0.93, bottom=0.09)
    # a: GC-matched exon - intron phyloP, same genes
    ax = fig.add_subplot(gs[0, 0])
    letter(ax, "a", x=-0.36)

    def gm(name, metric="phylop_p447"):
        t = tsv("comp", f"{name}_paired_summary.tsv").set_index(["metric", "comparison"]).loc[(metric, "gc_matched")]
        return t.exon_minus_intron, t.ci_lo, t.ci_hi
    rows = [("protein-coding", *gm("pc"), C["pc"], "s", True), ("lncRNA, ordinary", *gm("lnc"), C["ord"], "s", True)]
    for g in ("AluS", "AluJ", "MIR", "L2abc"):
        rows.append((f"lncRNA, {LABEL[g]}", *gm(f"te_{g}"), C[g], "o", True))
    ladder(ax, rows, (-0.1, 0.16), "exon − intron, mean phyloP (447 mammals)")
    ax.set_title("GC-matched conservation, same genes")

    # b: upper limits on TE-specific change in variant frequency (TE minus ordinary sequence, same genes)
    axb = fig.add_subplot(gs[0, 1])
    letter(axb, "b", x=-0.36)
    # gene-level excess (TE exon-intron minus ordinary exon-intron), variable callable non-CpG sites;
    # values from 12_gene_trend.py --baseline (work/trend/te_vs_nonte_genes.tsv), bootstrap over genes
    ex = {"AluS": {"pooled": (-0.0006, -0.0057, 0.0046), "sense": (-0.0041, -0.0118, 0.0045), "antisense": (0.0006, -0.0060, 0.0075)},
          "AluJ": {"pooled": (-0.0054, -0.0126, 0.0019), "sense": (0.0089, -0.0046, 0.0231), "antisense": (-0.0094, -0.0190, 0.0000)},
          "MIR": {"pooled": (0.0014, -0.0040, 0.0068), "sense": (-0.0009, -0.0093, 0.0074), "antisense": (0.0012, -0.0063, 0.0088)},
          "L2abc": {"pooled": (0.0004, -0.0063, 0.0069), "sense": (0.0031, -0.0065, 0.0124), "antisense": (-0.0006, -0.0102, 0.0088)}}
    base = {}
    for g in ex:
        b = tsv("comp", f"te_{g}_paired_bins.tsv")
        base[g] = (b.poly_c_noncpg_intron * b.units_intron).sum() / b.units_intron.sum()
    order = ["AluS", "AluJ", "MIR", "L2abc"]
    y = np.arange(len(order))[::-1]
    axb.axvline(0, color=GRID, lw=0.8, zorder=0)
    for yi, g in zip(y, order):
        e, lo, hi = [100 * v / base[g] for v in ex[g]["pooled"]]
        axb.plot([-12, lo], [yi, yi], color=C[g], lw=5, alpha=0.18, solid_capstyle="butt")
        axb.plot([lo, hi], [yi, yi], color=C[g], lw=1.6)
        axb.plot(e, yi, "o", color=C[g], ms=4.6, zorder=3)
        for off, ori, mk in ((0.25, "sense", "^"), (-0.25, "antisense", "v")):
            e2, lo2, hi2 = [100 * v / base[g] for v in ex[g][ori]]
            axb.plot([lo2, hi2], [yi + off] * 2, color=C[g], lw=0.7, alpha=0.6)
            axb.plot(e2, yi + off, mk, color=C[g], ms=2.8, alpha=0.8)
        axb.text(lo - 0.3, yi - 0.36, f"≥ {lo:.1f}%", ha="right", va="center", fontsize=5.8, color=C[g])
    axb.set_yticks(y)
    axb.set_yticklabels([LABEL[g] for g in order])
    axb.tick_params(axis="y", length=0)
    axb.set_xlim(-12, 8)
    axb.set_ylim(-0.7, len(order) + 0.35)
    axb.set_xlabel("TE-specific change in variant frequency (%)")
    axb.legend(handles=[Line2D([], [], marker="o", ls="-", color="#777777", ms=4, label="both orient."),
                        Line2D([], [], marker="^", ls="", color="#777777", ms=3, label="sense"),
                        Line2D([], [], marker="v", ls="", color="#777777", ms=3, label="antisense"),
                        Patch(color="#777777", alpha=0.18, label="excluded")],
               loc="upper center", ncol=4, fontsize=5.6, handlelength=1.1, columnspacing=0.7,
               handletextpad=0.3, borderaxespad=0.0)
    axb.set_title("Upper limits, TE minus ordinary sequence")

    # c: centered profiles (window minus the element's median window), same genes
    sub = GridSpecFromSubplotSpec(2, 2, subplot_spec=gs[1, :], width_ratios=[320, 280], hspace=0.7, wspace=0.1)
    layers = [("population", "variants", "RdBu_r", 0.25), ("phylop_p447", "phyloP", "RdBu", 0.35)]
    first = None
    for li, (layer, lname, cmap, vmax) in enumerate(layers):
        for gi, (g, run, L) in enumerate((("AluS", "AluS_b20_paired", 300), ("MIR", "MIR_b20_paired", 262))):
            axh = fig.add_subplot(sub[li, gi])
            if first is None:
                first = axh
                letter(axh, "c", x=-0.2, y=1.1)
            pr = tsv("prof_v5", run, "profiles.tsv")
            pr = pr[pr.layer == layer]
            nb = int(pr["bin"].max()) + 1
            M = np.full((2, nb), np.nan)
            Q = np.full((2, nb), np.nan)
            for ri, ori in enumerate(("sense", "antisense")):
                for r in pr[pr.orientation == ori].itertuples():
                    M[ri, int(r.bin)] = r.effect_centered
                    Q[ri, int(r.bin)] = r.q_centered
            im = axh.imshow(np.ma.masked_invalid(M), aspect="auto", cmap=cmap, vmin=-vmax, vmax=vmax,
                            extent=(0, nb * 20, 1.5, -0.5), interpolation="nearest")
            im.cmap.set_bad("#EEEEEE")
            for ri in range(2):
                for bi in range(nb):
                    if np.isfinite(Q[ri, bi]) and Q[ri, bi] < 0.1:
                        axh.plot(bi * 20 + 10, ri, "o", ms=3.2, mfc="none", mec=INK, mew=0.9)
            axh.set_yticks([0, 1])
            axh.set_yticklabels(["sense", "antisense"] if gi == 0 else ["", ""])
            axh.tick_params(axis="y", length=0)
            axh.set_xlim(0, nb * 20)
            for sp in axh.spines.values():
                sp.set_visible(False)
            if li == 0:
                axh.set_title(f"{LABEL[g]}, centered profile (same genes)", fontsize=6.8)
            if li == 1 and gi == 0:
                axh.set_xlabel("consensus position")
            if gi == 0:
                axh.text(-0.2, 0.5, lname, transform=axh.transAxes, rotation=90, va="center", ha="center",
                         fontsize=6.5, fontweight="bold")
            if gi == 1:
                cax = axh.inset_axes([1.03, 0.0, 0.025, 1.0])
                cb = fig.colorbar(im, cax=cax)
                cb.ax.tick_params(labelsize=5.5, length=1.5)
                cb.outline.set_linewidth(0.4)
                cb.set_label("vs median" if li == 0 else "vs median", fontsize=5.5)
    fig.text(0.9, 0.035, "window minus median window; blue = more constrained; circles: q < 0.1",
             ha="right", fontsize=5.8, color="#555555")
    save(fig, "Fig3_no_te_constraint")


# ----------------------------------------------------------------------------- Figure 4
def fig4():
    fig = plt.figure(figsize=(180 * MM, 74 * MM))
    gs = GridSpec(1, 3, figure=fig, width_ratios=[1, 1.05, 1.0], wspace=0.78, left=0.09, right=0.92, top=0.8,
                  bottom=0.22)
    # a: conserved stretches by subfamily and compartment
    ax = fig.add_subplot(gs[0, 0])
    letter(ax, "a", x=-0.3, y=1.22)
    cs = tsv("side", "conserved_summary.tsv")
    fams = ["MIR", "MIRb", "MIRc", "MIR3", "L2a", "L2b", "L2c", "L2d"]
    y = np.arange(len(fams))[::-1]
    for yi, f in zip(y, fams):
        t = cs[cs.family == f].set_index("compartment")
        for comp, col, off in (("exon", C["exon"], 0.22), ("intron", C["intron"], 0.0), ("intergenic", C["intergenic"], -0.22)):
            if comp in t.index:
                ax.barh(yi + off, 100 * t.loc[comp, "share_gt1.5"], height=0.22, color=col)
    ax.set_yticks(y)
    ax.set_yticklabels(fams)
    ax.tick_params(axis="y", length=0)
    ax.set_xlabel("copies with a conserved 20-nt stretch (%)")
    ax.legend(handles=[Patch(color=C["exon"], label="exon"), Patch(color=C["intron"], label="intron"),
                       Patch(color=C["intergenic"], label="intergenic")], loc="lower center",
              bbox_to_anchor=(0.45, 1.0), ncol=3, handlelength=1, columnspacing=0.8, borderaxespad=0.2)
    ax.set_title("Conserved fragments are not\nexon-specific", pad=14)

    # b: localization: exonic vs intronic TE density
    axb = fig.add_subplot(gs[0, 1])
    letter(axb, "b", x=-0.75, y=1.22)
    runs = [("both lines", "localization_models.tsv"), ("HepG2", "loc_HepG2_models.tsv"), ("K562", "loc_K562_models.tsv")]
    cats = [("L2b_MIRb_MIRc", "L2b+MIRb+MIRc", C["MIR"]), ("L2_all", "all L2", C["L2abc"]), ("MIR_all", "all MIR", "#7FBF7F"),
            ("Alu", "Alu", C["AluS"])]
    rows = []
    for cat, clab, col in cats:
        for rlab, f in runs:
            m = tsv("side", f)
            m = m[(m.category == cat) & (m.orientation == "all")].set_index("term")
            rows.append((f"{clab}, {rlab}", m.loc["exonic_density"], m.loc["intronic_density"], col))
    n = len(rows)
    y = np.arange(n)[::-1]
    axb.axvline(0, color=GRID, lw=0.8, zorder=0)
    for yi, (lab, E, I, col) in zip(y, rows):
        for r, off, mfc in ((E, 0.17, col), (I, -0.17, "white")):
            axb.plot([r.coef_per_sd - 1.96 * r.se, r.coef_per_sd + 1.96 * r.se], [yi + off] * 2, color=col, lw=1.0)
            axb.plot(r.coef_per_sd, yi + off, "o" if mfc != "white" else "s", ms=3.2, color=col, mfc=mfc, mew=0.9)
    for k in range(1, len(cats)):
        axb.axhline(n - 3 * k - 0.5, color="#EEEEEE", lw=0.6)
    axb.set_yticks(y)
    axb.set_yticklabels([r[0] for r in rows], fontsize=5.8)
    axb.tick_params(axis="y", length=0)
    axb.set_ylim(-0.7, n - 0.3)
    axb.set_xlabel("effect on log2 nucleus/cytoplasm\n(per SD of TE density)")
    axb.legend(handles=[Line2D([], [], marker="o", ls="", color="#555555", ms=3.2, label="exonic density"),
                        Line2D([], [], marker="s", ls="", mfc="white", mec="#555555", ms=3.2, label="intronic density")],
               loc="lower center", bbox_to_anchor=(0.45, 1.0), ncol=2, fontsize=5.8, borderaxespad=0.2,
               handletextpad=0.3)
    axb.set_title("Introns predict localization as\nwell as exons", pad=14)

    # c: RNA-property trends of the TE-specific excess
    axc = fig.add_subplot(gs[0, 2])
    letter(axc, "c", x=-0.75, y=1.22)
    tr = tsv("trend", "te_vs_nonte_trend.tsv")
    props = [("stab", "stability"), ("nuc", "nuclear"), ("splice", "splicing"), ("synth", "synthesis"),
             ("gtex", "expression"), ("tau", "tissue spec.")]
    grp = [("AluS*", "Alu S"), ("AluJ*", "Alu J"), ("MIR*", "MIR"), ("L2abc", "L2 a/b/c")]
    blocks = [("excess_phylop_p447", "phyloP"), ("excess_poly_noncpg", "variants")]
    M = np.full((len(blocks) * len(grp), len(props)), np.nan)
    P = np.full_like(M, np.nan)
    labels = []
    for bi, (layer, bl) in enumerate(blocks):
        for gi, (g, gl) in enumerate(grp):
            labels.append(f"{gl} · {bl}")
            for pi, (pp, _) in enumerate(props):
                r = tr[(tr.subfamily == g) & (tr.orientation == "pooled") & (tr.layer == layer) & (tr.property == pp)]
                if len(r):
                    M[bi * len(grp) + gi, pi] = r.spearman_rho.iloc[0]
                    P[bi * len(grp) + gi, pi] = r.p.iloc[0]
    im = axc.imshow(M, cmap="PuOr", vmin=-0.3, vmax=0.3, aspect="auto")
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            if np.isfinite(P[i, j]) and P[i, j] < 0.05:
                axc.plot(j, i, "o", ms=2.6, color=INK)
    axc.axhline(len(grp) - 0.5, color="white", lw=2)
    axc.set_xticks(range(len(props)))
    axc.set_xticklabels([p[1] for p in props], rotation=40, ha="right", fontsize=5.8)
    axc.set_yticks(range(len(labels)))
    axc.set_yticklabels(labels, fontsize=5.8)
    axc.tick_params(length=0)
    for sp in axc.spines.values():
        sp.set_visible(False)
    cax = axc.inset_axes([1.04, 0.0, 0.05, 1.0])
    cb = fig.colorbar(im, cax=cax)
    cb.ax.tick_params(labelsize=5.5, length=1.5)
    cb.outline.set_linewidth(0.4)
    cb.set_label("Spearman ρ across genes", fontsize=5.5)
    axc.set_title("No dose effect of RNA\nproperties (dots: P < 0.05)", pad=14)
    save(fig, "Fig4_fragments_rna")


# ----------------------------------------------------------------------------- Figure 5
def fig5():
    fig = plt.figure(figsize=(180 * MM, 100 * MM))
    gs = GridSpec(2, 2, figure=fig, width_ratios=[1.15, 1], height_ratios=[1, 1], wspace=0.42, hspace=0.62,
                  left=0.07, right=0.98, top=0.9, bottom=0.08)
    # a: stacked shares of features inside TEs by class
    ax = fig.add_subplot(gs[0, 0])
    letter(ax, "a", x=-0.14)
    af = tsv("side", "architecture_features.tsv").set_index(["gene_set", "feature"])
    feats = [("bg_exon", "exon\nsequence"), ("TSS", "start"), ("donor", "donor"), ("acceptor", "acceptor"), ("TES", "end")]
    classes = [("Alu", C["AluS"]), ("MIR", C["MIR"]), ("L2", C["L2abc"]), ("L1", "#8C564B"), ("LTR", "#56B4E9"),
               ("DNA", "#E7D43A")]
    x = np.arange(len(feats))
    for gi, (gset, off, lab) in enumerate((("lnc", -0.2, "lncRNA"), ("pc", 0.2, "protein-coding"))):
        bottom = np.zeros(len(feats))
        for cl, col in classes:
            v = np.array([af.loc[(gset, f), cl] for f, _ in feats]) * 100
            ax.bar(x + off, v, width=0.36, bottom=bottom, color=col, edgecolor="white", lw=0.3,
                   alpha=1.0 if gset == "lnc" else 0.55)
            bottom += v
        other = np.array([af.loc[(gset, f), "any TE"] for f, _ in feats]) * 100 - bottom
        ax.bar(x + off, other, width=0.36, bottom=bottom, color="#BBBBBB", edgecolor="white", lw=0.3,
               alpha=1.0 if gset == "lnc" else 0.55)
        tot = bottom + other
        for xi, t in zip(x + off, tot):
            ax.text(xi, t + 0.8, f"{t:.0f}", ha="center", fontsize=5.6)
    ax.set_xticks(x)
    ax.set_xticklabels([f[1] for f in feats], fontsize=6.3)
    ax.set_ylabel("inside a TE (%)")
    ax.set_ylim(0, 48)
    ax.set_xlabel("left bar: lncRNA · right (pale): protein-coding", fontsize=6, color="#555555")
    ax.legend(handles=[Patch(color=c, label=l) for l, c in classes] + [Patch(color="#BBBBBB", label="other")],
              loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=7, fontsize=5.6, handlelength=0.9,
              columnspacing=0.6, handletextpad=0.3, borderaxespad=0.2)
    ax.set_title("TEs supply lncRNA starts, splice sites and ends", pad=16)

    # b: LTR subfamilies at lncRNA starts
    axb = fig.add_subplot(gs[0, 1])
    letter(axb, "b", x=-0.3)
    lt = tsv("splice", "ltr_tss.tsv").sort_values("lnc_TSS_n", ascending=False).head(12)
    famcol = {"ERV1": "#56B4E9", "ERVL-MaLR": "#1F5F8B", "ERVL": "#8FB8DE", "ERVK": "#B3D4EA"}
    y = np.arange(len(lt))[::-1]
    for yi, r in zip(y, lt.itertuples()):
        axb.barh(yi, r.lnc_TSS_n, color=famcol.get(r.rfamily, "#999999"), height=0.7)
        axb.text(r.lnc_TSS_n + 8, yi, f"×{r.lnc_TSS_obs_exp:.1f} vs exons; ×{r.lnc_vs_pc_TSS:.0f} vs PC",
                 va="center", fontsize=5.4, color="#444444")
    axb.set_yticks(y)
    axb.set_yticklabels(lt.name)
    axb.tick_params(axis="y", length=0)
    axb.set_xlim(0, 760)
    axb.set_xlabel("lncRNA transcription starts")
    axb.legend(handles=[Patch(color=c, label=f) for f, c in famcol.items() if f in set(lt.rfamily)],
               loc="lower right", fontsize=5.8, handlelength=1, borderaxespad=0)
    axb.set_title("Most used LTRs at lncRNA starts")

    # c: Alu splice-site hotspots with positions reported for protein-coding Alu exons
    sub = GridSpecFromSubplotSpec(1, 2, subplot_spec=gs[1, :], wspace=0.18)
    hs = tsv("side", "architecture_hotspots.tsv")
    hs = hs[(hs.gene_set == "lnc") & (hs.rfamily == "Alu")]
    lit = {"antisense": {"acceptor": [116, 255, 275, 279], "donor": [23, 158]},
           "sense": {"acceptor": [], "donor": [51, 69, 101]}}
    for ci, ori in enumerate(("antisense", "sense")):
        axh = fig.add_subplot(sub[0, ci])
        if ci == 0:
            letter(axh, "c", x=-0.12, y=1.08)
        bins = np.arange(0, 316, 3)
        for feat, sign, col in (("acceptor", 1, "#C0392B"), ("donor", -1, "#2E6DA4")):
            d = hs[(hs.orientation == ori) & (hs.feature == feat)]
            h, _ = np.histogram(d.cons_pos, bins=bins, weights=d.n)
            axh.bar(bins[:-1] + 1.5, sign * h, width=3, color=col, lw=0)
            top = (h.max() if len(h) else 1) * 1.08
            for p in lit[ori][feat]:
                axh.plot(p, sign * top, marker="v" if sign > 0 else "^", ms=3.6, color=INK)
        axh.axhline(0, color=INK, lw=0.5)
        axh.set_xlim(0, 315)
        axh.set_xlabel("Alu consensus position")
        yl = axh.get_ylim()
        axh.set_yticks(axh.get_yticks())
        axh.set_yticklabels([f"{abs(int(t))}" for t in axh.get_yticks()])
        axh.set_ylim(yl)
        if ci == 0:
            axh.set_ylabel("lncRNA splice sites")
        axh.set_title(f"Alu copies {ori} to the transcript", fontsize=6.8)
        if ci == 1:
            axh.legend(handles=[Patch(color="#C0392B", label="acceptors"), Patch(color="#2E6DA4", label="donors"),
                                Line2D([], [], marker="v", ls="", color=INK, ms=3.6,
                                       label="sites of Alu exons in\nprotein-coding genes")],
                       loc="lower right", fontsize=5.8, handlelength=1, borderaxespad=0.2)
    save(fig, "Fig5_architecture")


# ----------------------------------------------------------------------------- Figure 6
def fig6():
    fig = plt.figure(figsize=(180 * MM, 70 * MM))
    gs = GridSpec(1, 3, figure=fig, wspace=0.2, left=0.23, right=0.985, top=0.86, bottom=0.2)
    st = tsv("splice", "res_stats.tsv")
    pick = [("protein-coding splice sites vs intronic flank", "all", "protein-coding (vs intron)", C["pc"]),
            ("lncRNA non-TE splice sites vs intronic flank", "all", "lncRNA, outside TEs (vs intron)", C["ord"])]
    te = "lncRNA TE-derived splice sites vs same positions in unused copies"
    tecol = "#C98A00"
    for sub, lab, col in (("all", "lncRNA, in TEs (vs unused copies)", tecol), ("Alu", "· Alu", C["AluS"]),
                          ("MIR", "· MIR", C["MIR"]), ("L2", "· L2", C["L2abc"]), ("donor", "· donors", tecol),
                          ("acceptor", "· acceptors", tecol),
                          ("antisense-Alu acceptor hotspots", "· antisense-Alu acceptor hotspots", C["AluY"])):
        pick.append((te, sub, lab, col))
    rows_v, rows_p, rows_s = [], [], []
    for comp, sub, lab, col in pick:
        r = st[(st.comparison == comp) & (st.subset == sub)].iloc[0]
        rows_v.append((lab, r.mh_log_or, r.mh_log_or_lo, r.mh_log_or_hi, col, "o", True))
        rows_s.append((lab, r.singleton_diff, r.singleton_diff_lo, r.singleton_diff_hi, col, "o", True))
        rows_p.append((lab, r.phylop_p447_diff, r.phylop_p447_diff_lo, r.phylop_p447_diff_hi, col, "o", True))
    for i, (rows, xl, xlab, title) in enumerate((
            (rows_v, (-0.5, 0.3), "log OR of a variable site\n(below 0 = fewer variants)", "Human variants"),
            (rows_s, (-0.1, 0.15), "difference in singleton share\n(above 0 = rarer variants)", "Rare variants"),
            (rows_p, (-0.25, 0.45), "difference in phyloP\n(above 0 = more conserved)", "Conservation (mammals)"))):
        ax = fig.add_subplot(gs[0, i])
        letter(ax, "abc"[i], x=-0.13 if i else -0.66, y=1.02)
        ladder(ax, rows, xl, xlab)
        ax.axhline(len(rows) - 2.5, color="#E5E5E5", lw=0.6)
        if i:
            ax.set_yticklabels([])
        ax.set_title(title)
    save(fig, "Fig6_splice_sites")


if __name__ == "__main__":
    from matplotlib.backends.backend_pdf import PdfPages
    want = set(sys.argv[3:]) if len(sys.argv) > 3 else None
    if want is None:
        PDF_ALL = PdfPages(os.path.join(OUT, "RIDL2_main_figures.pdf"))
    for k, f in (("1", fig1), ("2", fig2), ("3", fig3), ("4", fig4), ("5", fig5), ("6", fig6)):
        if want is None or k in want:
            f()
    if PDF_ALL is not None:
        PDF_ALL.close()
