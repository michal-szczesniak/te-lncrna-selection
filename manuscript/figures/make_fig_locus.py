#!/usr/bin/env python3
"""Example loci: lncRNA exons whose splice acceptor comes from an antisense Alu copy.

Reads results/locus/bundle.json of step 44 (patch 9) and draws, for each chosen site, the gene model with the
Alu copy and the GTEx junction reads, and the sequence around the acceptor in human (transcript orientation, as
RNA) aligned to the Alu consensus and to the orthologous sequence of five primates and mouse.

    python make_fig_locus.py <bundle.json> <out_dir> <site_id> [<site_id> ...]

Saves <out_dir>/Fig5.png and .pdf, or <out_dir>/$FIGNAME.png and .pdf when FIGNAME is set.
"""
import os
import sys
import glob

import json

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
from matplotlib.gridspec import GridSpec
from matplotlib.patches import FancyArrow, Rectangle

for _f in glob.glob(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "fonts", "*.ttf")):
    fm.fontManager.addfont(_f)
plt.rcParams.update({"font.family": "Liberation Sans", "font.size": 7.5, "axes.titlesize": 8,
                     "axes.titleweight": "bold", "axes.titlelocation": "left", "pdf.fonttype": 42})
MM = 1 / 25.4
INK, GREY = "#222222", "#8C8C8C"
EXON, ALU, ACC, LOST = "#4D4D4D", "#E69F00", "#D55E00", "#0072B2"
ACC_TXT, ALU_TXT = "#A33A00", "#8A5A00"  # darker shades for text on white
SPECIES = {"panTro6": "chimpanzee", "gorGor6": "gorilla", "ponAbe3": "orangutan", "rheMac10": "macaque",
           "calJac4": "marmoset", "mm39": "mouse"}
UP, DOWN = 24, 10  # nt shown upstream of the AG and downstream of it (exon)


def letter(ax, s, x=-0.06, y=1.02):
    ax.text(x, y, s, transform=ax.transAxes, fontsize=10, fontweight="bold", va="bottom", ha="left", color=INK)


def gene_panel(ax, b, title):
    s, cp = b["site"], b["copy"]
    strand = s["tx_strand"]
    tx = sorted(b["transcripts"], key=lambda t: (not t["canonical"], t["transcript_id"]))
    ik = s["site_id"][:-1]
    i0, i1 = [int(x) for x in ik.split(":")[1][:-1].split("-")]
    acc_pos = i1 if strand == "+" else i0  # exon edge made by the acceptor
    donor = i0 if strand == "+" else i1    # the upstream donor that splices to the Alu acceptor
    # one transcript: the canonical one if it uses the Alu acceptor, otherwise the first one that does
    use = [t for t in tx if any(e[0] == acc_pos or e[1] == acc_pos for e in t["exons"])]
    t = next((u for u in use if u["canonical"]), use[0])
    ex = t["exons"]
    lo, hi = min(ex[0][0], cp["start"]), max(ex[-1][1], cp["end"])
    pad = 0.03 * (hi - lo)
    ax.set_xlim(lo - pad, hi + pad)
    ax.plot([ex[0][0], ex[-1][1]], [0, 0], color=GREY, lw=0.8, zorder=1)
    for a, e in ex:
        ax.add_patch(Rectangle((a, -0.22), e - a, 0.44, facecolor=EXON, edgecolor=EXON, lw=0.4, zorder=2))
    # GTEx junctions, sashimi style: introns of this transcript above the line, other junctions between its
    # exons below; width and shade by reads, scaled to the most used junction of this gene
    introns = [(ex[k][1], ex[k + 1][0]) for k in range(len(ex) - 1)]
    J = {(j["start"], j["end"]): j for j in b["junctions"]}
    edges = {e for a_e in ex for e in a_e}
    above = [J[k] for k in introns if k in J and (J[k].get("total_reads") or 0) > 0]
    below = [j for k, j in J.items() if k not in introns and k[0] in edges and k[1] in edges
             and (j.get("total_reads") or 0) > 0]
    rmax = max([1.0] + [j["total_reads"] for j in above + below])
    top = max([1.0] + [np.log10(j["total_reads"]) for j in above + below])
    for side, js in ((1, above), (-1, below)):
        for j in js:
            x0, x1, r = j["start"], j["end"], j["total_reads"]
            f = r / rmax
            h = 0.35 + (0.9 if side > 0 else 0.5) * np.log10(r) / top
            xs = np.linspace(x0, x1, 60)
            ys = side * (0.25 + h * np.sin(np.linspace(0, np.pi, 60)))
            alu_j = (x0, x1) == (i0, i1)
            ax.plot(xs, ys, color=ACC if alu_j else str(round(0.72 - 0.4 * f, 3)), lw=0.5 + 2.1 * f, zorder=1,
                    solid_capstyle="butt")
            yl = ys.max() + 0.05 + 0.06 * f if side > 0 else ys.min() - 0.05 - 0.06 * f
            ax.text((x0 + x1) / 2, yl, f"{int(r):,}", ha="center", va="bottom" if side > 0 else "top",
                    fontsize=5.4, color=ACC_TXT if alu_j else "#444444", fontweight="bold" if alu_j else "normal",
                    zorder=3, bbox=dict(facecolor="white", edgecolor="none", pad=0.4))
    for k in introns:  # introns of this transcript that the GTEx junction table does not contain
        if k not in J:
            ax.text((k[0] + k[1]) / 2, 0.3, "n/a", ha="center", va="bottom", fontsize=5.0, color=GREY,
                    style="italic")
    # Alu copy, drawn in its own orientation (antisense to the transcript)
    ya = -2.05 if below else -0.75
    L = cp["end"] - cp["start"]
    hl = min(L * 0.35, 0.02 * (hi - lo))
    if cp["te_strand"] == "+":
        ax.add_patch(FancyArrow(cp["start"], ya, L, 0, width=0.32, head_width=0.5, head_length=hl,
                                length_includes_head=True, color=ALU, lw=0))
    else:
        ax.add_patch(FancyArrow(cp["end"], ya, -L, 0, width=0.32, head_width=0.5, head_length=hl,
                                length_includes_head=True, color=ALU, lw=0))
    ax.text((cp["start"] + cp["end"]) / 2, ya - 0.45, f"{cp['subfamily_orig']} (antisense)", ha="center",
            va="top", fontsize=6.2, color=ALU_TXT)
    ax.plot([acc_pos, acc_pos], [ya + 0.25, 1.9], color=ACC, lw=0.8, ls=(0, (2, 2)), zorder=0)
    reads, n5 = int(s["total_reads"]), int(s["samples_ge5"])
    dshare = J[(i0, i1)].get("donor_share") if (i0, i1) in J else None
    skip = sum(j.get("donor_share") or 0 for j in below if donor in (j["start"], j["end"]))
    tissue = str(s.get("max_tissue") or "").split(" - ")[0].lower()
    ax.text(lo - pad, ya - 1.15, f"Alu acceptor (dashed): {reads:,} GTEx reads, \u22655 reads in {n5:,} samples",
            ha="left", va="top", fontsize=6.0, color=INK)
    use_txt = "Reads from the upstream donor: " + (f"{100 * dshare:.0f}% here" if dshare is not None
                                                     and np.isfinite(dshare) else "n/a")
    if skip > 0:
        use_txt += f", {100 * skip:.0f}% skip the exon (lower arc)"
    ax.text(lo - pad, ya - 1.6, use_txt, ha="left", va="top", fontsize=6.0, color=INK)
    ax.text(lo - pad, ya - 2.05, f"Most reads per sample in {tissue}", ha="left", va="top", fontsize=6.0, color=INK)
    span = hi - lo
    sb = 10 ** np.floor(np.log10(span * 0.25))
    sb = sb * (5 if span * 0.25 >= 5 * sb else 2 if span * 0.25 >= 2 * sb else 1)
    ax.plot([hi - sb, hi], [ya - 2.65, ya - 2.65], color=INK, lw=0.8)
    ax.text(hi - sb / 2, ya - 2.75, f"{sb / 1000:g} kb", ha="center", va="top", fontsize=6)
    ax.text(lo - pad, 2.35, f"{b['gene_name']}  \u00b7  {s['chrom']}, {strand} strand  \u00b7  "
            f"{t['transcript_id'].split('.')[0]}" + (" (canonical)" if t["canonical"] else ""),
            ha="left", va="bottom", fontsize=6.8, fontweight="bold")
    ax.set_ylim(-5.15, 2.8)  # same vertical scale in every panel
    ax.axis("off")
    ax.set_title("     " + title, pad=6)


def seq_panel(ax, b):
    w = b["window"]
    k0 = len(w["human"]) // 2 - 1  # first base of the AG (flank nt on each side)
    sl = slice(k0 - UP, k0 + 2 + DOWN)
    rows = [("human", w["human"][sl], INK)]
    rows.append(("Alu consensus", "".join(c if c != " " else "." for c in w["consensus"][sl]), ALU))
    for sp in b["species_order"]:
        rows.append((SPECIES.get(sp, sp), w["species"].get(sp, "")[sl], INK))
    n = UP + 2 + DOWN
    ax.set_xlim(-0.5, n + 7.5)
    ax.set_ylim(len(rows) - 0.4, -1.2)
    ax.axis("off")
    hum = rows[0][1]
    for i, (lab, seq, col) in enumerate(rows):
        ax.text(-0.9, i, lab, ha="right", va="center", fontsize=6.6, color=col)
        if set(seq) <= {"-"}:
            ax.text(n / 2, i, "not aligned (no Alu copy at this position)", ha="center", va="center",
                    fontsize=6.2, color=GREY, style="italic")
            continue
        for j, c in enumerate(seq):
            r = "U" if c == "T" else c
            same = i == 0 or c == hum[j] or lab == "Alu consensus"
            colr = col if same else LOST
            fw = "bold" if (UP <= j < UP + 2) else "normal"
            ax.text(j, i, r, ha="center", va="center", fontsize=6.6, family="DejaVu Sans Mono", color=colr,
                    fontweight=fw)
        if i >= 2:
            ag = seq[UP:UP + 2]
            mark = "AG" if ag == "AG" else ("no AG" if "-" not in ag else "")
            ax.text(n + 0.8, i, mark, ha="left", va="center", fontsize=6.2,
                    color=INK if ag == "AG" else LOST)
    ax.add_patch(Rectangle((UP - 0.5, -0.5), 2, len(rows), facecolor="#FCE6CF", edgecolor="none", zorder=0))
    ax.plot([UP + 1.5] * 2, [-0.9, len(rows) - 0.5], color=ACC, lw=0.8, ls=(0, (2, 2)))
    ax.text(UP / 2, -1.0, "intron (U-rich tract)", ha="center", va="center", fontsize=6.2, color="#555555")
    ax.text(UP + 2 + DOWN / 2, -1.0, "exon", ha="center", va="center", fontsize=6.2, color="#555555")


def main():
    bundle, out = sys.argv[1], sys.argv[2]
    ids = sys.argv[3:]
    B = {b["site"]["site_id"]: b for b in json.load(open(bundle))}
    os.makedirs(out, exist_ok=True)
    fig = plt.figure(figsize=(180 * MM, 62 * MM * len(ids)))
    gs = GridSpec(len(ids), 2, figure=fig, width_ratios=[1.05, 1], wspace=0.28, hspace=0.45,
                  left=0.03, right=0.99, top=0.95, bottom=0.04)
    titles = {"ape_insertion": "An Alu shared only with the great apes provides the acceptor",
              "lost_in_monkey": "The Alu is older, but macaque lacks the acceptor AG"}
    for i, sid in enumerate(ids):
        b = B[sid]
        ax = fig.add_subplot(gs[i, 0])
        letter(ax, "AB"[i], x=-0.03, y=1.06)
        gene_panel(ax, b, titles.get(b["site"]["locus_class"], ""))
        ax2 = fig.add_subplot(gs[i, 1])
        seq_panel(ax2, b)
    name = os.environ.get("FIGNAME", "Fig5")  # e.g. FIGNAME=Fig5_v3 for a test render with its own file name
    fig.savefig(os.path.join(out, name + ".png"), dpi=300, facecolor="white")
    fig.savefig(os.path.join(out, name + ".pdf"), facecolor="white")


if __name__ == "__main__":
    main()
