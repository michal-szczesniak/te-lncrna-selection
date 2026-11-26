#!/usr/bin/env python3
"""Supplementary Fig. S8: lncRNA transcription starts in LTR elements (step 42).
(A) Share of starts in LTRs that lie in the same orientation as the transcript, per LTR group, with a 95% binomial
interval, against the share expected from the orientation of LTR sequence in the exons of the same genes.
(B) Starts in sense LTRs along the consensus (10-nt bins) for each named LTR group; shaded, the most used 50-nt
window; text, its share of the group's sense starts and the share expected for the most used window if starts were
placed at random in the same copies. The x axis of each panel ends at the 98th percentile of positions.
Usage: python make_figS8.py <results/revision6> <out stem>"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
import numpy as np
import pandas as pd
from scipy.stats import beta

# Liberation Sans ships with the code (manuscript/fonts, SIL OFL 1.1), so figures look the same on any machine
import glob as _glob
from matplotlib import font_manager as _fm
for _f in _glob.glob(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "fonts", "*.ttf")):
    _fm.fontManager.addfont(_f)
D, OUT = sys.argv[1], sys.argv[2]
plt.rcParams.update({
    "font.family": "Liberation Sans", "font.size": 7.5, "axes.titlesize": 8, "axes.labelsize": 7.5,
    "xtick.labelsize": 6.5, "ytick.labelsize": 6.5, "legend.fontsize": 6.8, "axes.linewidth": 0.6,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6, "xtick.major.size": 2.5, "ytick.major.size": 2.5,
    "axes.spines.top": False, "axes.spines.right": False, "pdf.fonttype": 42, "axes.titleweight": "bold",
    "axes.titlelocation": "left", "legend.frameon": False})
INK, GRID, SKY, VERM, SHADE = "#222222", "#d9d9d9", "#56B4E9", "#D55E00", "#E6E6E6"
MM = 1 / 25.4


def letter(fig, ax, s, dx=-0.02):
    bb = ax.get_position()
    fig.text(bb.x0 + dx, bb.y1 + 0.035, s, fontsize=10, fontweight="bold", va="bottom", ha="right", color=INK)


S = pd.read_csv(os.path.join(D, "tss_ltr_summary.tsv"), sep="\t")
H = pd.read_csv(os.path.join(D, "tss_ltr_hist.tsv"), sep="\t")
S = S[S.starts >= 20]
named = S[~S.group.isin(["all LTR", "other LTR"])].sort_values("starts", ascending=False)
order = pd.concat([S[S.group == "all LTR"], named, S[S.group == "other LTR"]])

fig = plt.figure(figsize=(180 * MM, 105 * MM))
gs = GridSpec(1, 2, figure=fig, width_ratios=[0.85, 2.0], wspace=0.22, left=0.13, right=0.985, top=0.86,
              bottom=0.2)

# ---- A: orientation
ax = fig.add_subplot(gs[0, 0])
y = np.arange(len(order))[::-1]
for yi, (_, r) in zip(y, order.iterrows()):
    lo = beta.ppf(0.025, r.sense, r.starts - r.sense + 1) if r.sense > 0 else 0.0
    hi = beta.ppf(0.975, r.sense + 1, r.starts - r.sense) if r.sense < r.starts else 1.0
    col = INK if r.group == "all LTR" else SKY
    ax.plot([100 * lo, 100 * hi], [yi, yi], color=col, lw=1.2, solid_capstyle="butt")
    ax.plot(100 * r.share_sense, yi, "o", color=col, ms=4)
    if np.isfinite(r.expected_share_sense):
        ax.plot([100 * r.expected_share_sense] * 2, [yi - 0.32, yi + 0.32], color=VERM, lw=1.4)
ax.set_yticks(y)
ax.set_yticklabels([f"{g} ({int(n):,})" for g, n in zip(order.group, order.starts)])
ax.tick_params(axis="y", length=0)
ax.axvline(50, color=GRID, lw=0.8, zorder=0)
ax.set_xlim(40, 100)
ax.set_ylim(-0.7, len(order) - 0.3)
ax.set_xlabel("starts in sense LTRs (%)")
ax.set_title("Orientation", pad=6)
ax.plot([], [], "o", color=SKY, ms=4, label="observed, 95% interval")
ax.plot([], [], "|", color=VERM, ms=8, mew=1.4, label="expected from exonic\nLTR sequence")
ax.legend(loc="upper left", bbox_to_anchor=(-0.55, -0.13), ncol=1, borderaxespad=0.2, handletextpad=0.4)

# ---- B: position along the consensus, named groups
n = len(named)
ncol = 3
nrow = int(np.ceil(n / ncol))
sub = GridSpecFromSubplotSpec(nrow, ncol, subplot_spec=gs[0, 1], hspace=0.75, wspace=0.32)
axes = []
for i, (_, r) in enumerate(named.iterrows()):
    a = fig.add_subplot(sub[i // ncol, i % ncol])
    axes.append(a)
    h = H[(H.group == r.group) & (H.orientation == "sense")].sort_values("cons_bin")
    if len(h):
        cum = h.starts.cumsum() / h.starts.sum()
        xmax = float(h.cons_bin[cum >= 0.98].min()) + 10
        if np.isfinite(r.window_end):
            xmax = max(xmax, float(r.window_end) + 30)
        xmax = 50 * np.ceil(xmax / 50)
        h = h[h.cons_bin < xmax]
        a.bar(h.cons_bin + 5, h.starts, width=10, color=SKY, lw=0)
        a.set_xlim(0, xmax)
    if np.isfinite(r.window_start):
        a.axvspan(r.window_start, r.window_end, color=SHADE, lw=0, zorder=0)
        a.text(0.98, 0.97, f"{100 * r.window_share:.0f}% in {int(r.window_start)}–{int(r.window_end)}\n"
               f"random: {100 * r.window_share_expected:.0f}%", transform=a.transAxes, ha="right", va="top",
               fontsize=6.2, color=INK, linespacing=1.15)
    a.set_title(f"{r.group} ({int(r.sense):,} sense starts)", fontsize=7, pad=3)
    a.set_ylim(0, a.get_ylim()[1] * 1.35)
    a.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(4, integer=True))
    a.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(3, integer=True))
    if i % ncol == 0:
        a.set_ylabel("starts")
    if i // ncol == nrow - 1 or i + ncol >= n:
        a.set_xlabel("consensus position")
fig.canvas.draw()
letter(fig, ax, "A", dx=-0.1)
letter(fig, axes[0], "B", dx=-0.035)
fig.text(axes[0].get_position().x0, ax.get_position().y1 + 0.075, "Position of starts in sense LTRs",
         fontsize=8, fontweight="bold", ha="left", va="bottom", color=INK)
for ext in ("png", "pdf"):
    fig.savefig(f"{OUT}.{ext}", dpi=300)
print("written", OUT)
