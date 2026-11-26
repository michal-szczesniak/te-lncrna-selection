#!/usr/bin/env python3
"""Supplemental Fig. S3B: planted constraint scanned along the element with genes resampled (step 33).
Usage: python make_figS3b.py <power_scan dir with <group>_<orientation>_power.tsv> <out png/pdf stem>"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

D, OUT = sys.argv[1], sys.argv[2]
# Liberation Sans ships with the code (manuscript/fonts, SIL OFL 1.1), so figures look the same on any machine
import glob as _glob
from matplotlib import font_manager as _fm
for _f in _glob.glob(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "fonts", "*.ttf")):
    _fm.fontManager.addfont(_f)
plt.rcParams.update({
    "font.family": "Liberation Sans", "font.size": 7.5, "axes.titlesize": 8, "axes.labelsize": 7.5,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 6.8, "axes.linewidth": 0.6,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6, "xtick.major.size": 2.5, "ytick.major.size": 2.5,
    "axes.spines.top": False, "axes.spines.right": False, "pdf.fonttype": 42, "axes.titleweight": "bold",
    "axes.titlelocation": "left", "legend.frameon": False})
INK = "#222222"
COL = {"AluS": "#E69F00", "AluJ": "#E69F00", "MIR": "#009E73", "L2abc": "#CC79A7"}
LAB = {"AluS": "AluS", "AluJ": "AluJ", "MIR": "MIR", "L2abc": "L2a/b/c"}
MM = 1 / 25.4


def shade(c, f):  # f = 0 -> white, 1 -> colour
    r = np.array(to_rgb(c))
    return tuple(1 - f * (1 - r))


groups = ["AluS", "AluJ", "MIR", "L2abc"]
fig, axes = plt.subplots(4, 2, figsize=(180 * MM, 175 * MM), sharey=True)
fig.subplots_adjust(left=0.1, right=0.98, top=0.93, bottom=0.07, hspace=0.55, wspace=0.12)
for r, g in enumerate(groups):
    for c, ori in enumerate(("sense", "antisense")):
        ax = axes[r, c]
        f = os.path.join(D, f"{g}_{ori}_power.tsv")
        t = pd.read_csv(f, sep="\t")
        reps = int(t.reps.max())
        for red, f_ in ((0.1, 0.45), (0.25, 0.72), (0.5, 1.0)):
            for phi, ls in ((1.0, "-"), (0.5, "--")):
                q = t[(t.reduction == red) & (t.fraction_units == phi)].sort_values("seg_start")
                ax.plot(q.seg_start + 9.5, q.power_window, ls, color=shade(COL[g], f_), lw=1.1)
        q = t[t.reduction == 0].sort_values("seg_start")
        ax.plot(q.seg_start + 9.5, q.power_window, ":", color="#777777", lw=0.9)
        ax.set_ylim(-0.03, 1.03)
        ax.set_title(f"{LAB[g]} {ori}", fontsize=7.5)
        if r == 3:
            ax.set_xlabel("centre of the planted 20-nt segment (consensus position)")
        if c == 0:
            ax.set_ylabel("detection rate")
        ax.axhline(0.8, color="#cccccc", lw=0.6, zorder=0)
h = [Line2D([], [], color=shade("#555555", f_), lw=1.3, label=f"−{int(100 * red)}% of variants")
     for red, f_ in ((0.1, 0.45), (0.25, 0.72), (0.5, 1.0))]
h += [Line2D([], [], color="#555555", ls="-", lw=1.1, label="in all exonic copies"),
      Line2D([], [], color="#555555", ls="--", lw=1.1, label="in half of them"),
      Line2D([], [], color="#777777", ls=":", lw=0.9, label="nothing planted")]
fig.legend(handles=h, loc="upper center", ncol=6, bbox_to_anchor=(0.54, 0.995), handlelength=2.0, columnspacing=1.2)
fig.text(0.015, 0.985, "B", fontsize=10, fontweight="bold", va="top")
fig.savefig(OUT + ".png", dpi=300, facecolor="white")
fig.savefig(OUT + ".pdf", facecolor="white")
print("reps", reps)
