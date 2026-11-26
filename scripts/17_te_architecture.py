#!/usr/bin/env python3
"""TEs in lncRNA gene architecture: transcription starts, splice sites and transcript ends inside TEs.

For each gene set (--gene-set name=lnc_genes.tsv.gz from step 1; kept genes = empty drop_reason), all
annotated transcripts of the kept genes are taken from the GTF and four kinds of features are collected
(unique positions per gene): TSS (first base), TES (last base), donor (last exonic base before an
intron) and acceptor (first exonic base after an intron). Each position is checked against RepeatMasker
.out (all records). Expectations come from random positions in the exons of the same genes (bg_exon,
length-weighted) and in their gene bodies (bg_genebody).

For features inside Alu, MIR and L2 elements the consensus position is computed by linear interpolation
of the RepeatMasker coordinates, separately for TEs in the same (sense) or opposite (antisense)
orientation as the transcript: this shows which positions of the elements are used as splice sites,
starts or ends ("hotspots").

Outputs: <out>_features.tsv (share of features inside each TE class, obs/exp vs bg_exon),
         <out>_hotspots.tsv (counts by family, feature, orientation, consensus position),
         <out>_hotspots.png
"""
import argparse

import numpy as np
import pandas as pd

from ridl_utils import log, open_text, read_gtf

CLASSES = {"any TE": lambda c, f: c in ("SINE", "LINE", "LTR", "DNA", "Retroposon"),
           "Alu": lambda c, f: f == "Alu", "MIR": lambda c, f: f == "MIR",
           "L1": lambda c, f: f == "L1", "L2": lambda c, f: f == "L2",
           "LTR": lambda c, f: c == "LTR", "DNA": lambda c, f: c == "DNA"}
HOT = ("Alu", "MIR", "L2")


def num(x):
    return int(x.strip("()"))


def read_rmsk(path, chroms):
    rows = []
    with open_text(path) as fh:
        for line in fh:
            f = line.split()
            if len(f) < 15 or not f[0].isdigit() or f[4] not in chroms:
                continue
            cf = f[10].split("/")
            plus = f[8] == "+"
            if plus:
                cs, ce = num(f[11]), num(f[12])
            else:
                ce, cs = num(f[12]), num(f[13])
            rows.append((f[4], int(f[5]) - 1, int(f[6]), "+" if plus else "-", cf[0], cf[-1], cs, ce))
    r = pd.DataFrame(rows, columns=["chrom", "start", "end", "strand", "rclass", "rfamily", "cs", "ce"])
    return {c: g.sort_values("start").reset_index(drop=True) for c, g in r.groupby("chrom")}


def annotate(pts, rm):
    """pts: DataFrame chrom, pos, strand (transcript strand). Adds rclass, rfamily, orientation, cons_pos."""
    out = []
    for c, g in pts.groupby("chrom"):
        g = g.copy()
        g["rclass"], g["rfamily"], g["orientation"], g["cons_pos"] = "", "", "", np.nan
        if c not in rm:
            out.append(g)
            continue
        r = rm[c]
        S, E = r["start"].to_numpy(), r["end"].to_numpy()
        P = g["pos"].to_numpy()
        idx = np.searchsorted(S, P, side="right") - 1
        hit = np.full(len(P), -1)
        for k in range(6):
            cand = idx - k
            ok = (hit < 0) & (cand >= 0)
            ok[ok] &= E[cand[ok]] > P[ok]
            hit[ok] = cand[ok]
        h = hit >= 0
        rr = r.iloc[hit[h]]
        g.loc[h, "rclass"] = rr["rclass"].to_numpy()
        g.loc[h, "rfamily"] = rr["rfamily"].to_numpy()
        g.loc[h, "orientation"] = np.where(rr["strand"].to_numpy() == g.loc[h, "strand"].to_numpy(),
                                           "sense", "antisense")
        frac = (P[h] - rr["start"].to_numpy()) / np.maximum(rr["end"].to_numpy() - rr["start"].to_numpy() - 1, 1)
        cs, ce = rr["cs"].to_numpy(), rr["ce"].to_numpy()
        g.loc[h, "cons_pos"] = np.round(np.where(rr["strand"].to_numpy() == "+", cs + frac * (ce - cs),
                                                 ce - frac * (ce - cs)))
        out.append(g)
    return pd.concat(out, ignore_index=True)


def features(tx, ex):
    ex = ex.merge(tx[["transcript_id", "strand"]], on="transcript_id", suffixes=("", "_t"))
    rows = []
    t = tx
    plus = t["strand"] == "+"
    rows.append(pd.DataFrame({"gene_id": t["gene_id"], "chrom": t["chrom"], "strand": t["strand"], "feature": "TSS",
                              "pos": np.where(plus, t["start"], t["end"] - 1)}))
    rows.append(pd.DataFrame({"gene_id": t["gene_id"], "chrom": t["chrom"], "strand": t["strand"], "feature": "TES",
                              "pos": np.where(plus, t["end"] - 1, t["start"])}))
    e = ex.sort_values(["transcript_id", "start"])
    nxt_start = e.groupby("transcript_id")["start"].shift(-1)
    m = nxt_start.notna()
    left = e[m]
    right_start = nxt_start[m].astype(np.int64)
    pl = left["strand"] == "+"
    # intron between left exon (ends at left.end) and next exon (starts at right_start)
    donor = np.where(pl, left["end"] - 1, right_start)
    acceptor = np.where(pl, right_start, left["end"] - 1)
    rows.append(pd.DataFrame({"gene_id": left["gene_id"], "chrom": left["chrom"], "strand": left["strand"],
                              "feature": "donor", "pos": donor}))
    rows.append(pd.DataFrame({"gene_id": left["gene_id"], "chrom": left["chrom"], "strand": left["strand"],
                              "feature": "acceptor", "pos": acceptor}))
    return pd.concat(rows, ignore_index=True).drop_duplicates(["gene_id", "feature", "pos"])


def backgrounds(tx, ex, genes, n, rng):
    rows = []
    em = ex.drop_duplicates(["gene_id", "start", "end"])
    for gid, g in em.groupby("gene_id"):
        L = (g["end"] - g["start"]).to_numpy()
        k = rng.choice(len(g), n, p=L / L.sum())
        pos = g["start"].to_numpy()[k] + (rng.random(n) * L[k]).astype(np.int64)
        rows.append(pd.DataFrame({"gene_id": gid, "chrom": g["chrom"].iloc[0], "strand": g["strand"].iloc[0],
                                  "feature": "bg_exon", "pos": pos}))
    gb = genes.set_index("gene_id")
    for gid in em["gene_id"].unique():
        s, e = int(gb.at[gid, "start"]), int(gb.at[gid, "end"])
        rows.append(pd.DataFrame({"gene_id": gid, "chrom": gb.at[gid, "chrom"], "strand": gb.at[gid, "strand"],
                                  "feature": "bg_genebody", "pos": rng.integers(s, e, n)}))
    return pd.concat(rows, ignore_index=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--gene-set", action="append", required=True, help="name=lnc_genes.tsv.gz (step 1)")
    ap.add_argument("--rmsk-out", required=True)
    ap.add_argument("--bg-per-gene", type=int, default=50)
    ap.add_argument("--chroms", default=",".join(f"chr{i}" for i in range(1, 23)))
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)
    chroms = a.chroms.split(",")
    genes, tx, ex = read_gtf(a.gtf, chroms)
    rm = read_rmsk(a.rmsk_out, set(chroms))
    log(f"repeat records: {sum(len(v) for v in rm.values())}")

    summ, hots = [], []
    for item in a.gene_set:
        name, path = item.split("=", 1)
        g = pd.read_csv(path, sep="\t")
        kept = set(g.loc[g["drop_reason"].fillna("") == "", "gene_id"])
        t = tx[tx["gene_id"].isin(kept)]
        e = ex[ex["gene_id"].isin(kept)]
        f = features(t, e)
        b = backgrounds(t, e, genes[genes["gene_id"].isin(kept)], a.bg_per_gene, rng)
        pts = annotate(pd.concat([f, b], ignore_index=True), rm)
        log(f"{name}: genes {len(kept)}, transcripts {len(t)}, features {len(f)}")
        for feat, s in pts.groupby("feature"):
            r = {"gene_set": name, "feature": feat, "n": len(s)}
            for cl, fn in CLASSES.items():
                r[cl] = np.mean([fn(c, fm) for c, fm in zip(s["rclass"], s["rfamily"])])
            summ.append(r)
        h = pts[pts["rfamily"].isin(HOT) & pts["feature"].isin(["TSS", "TES", "donor", "acceptor"])]
        hots.append(h.groupby(["rfamily", "feature", "orientation", "cons_pos"]).size()
                    .rename("n").reset_index().assign(gene_set=name))
    s = pd.DataFrame(summ)
    for name, g in s.groupby("gene_set"):
        bg = g[g["feature"] == "bg_exon"].iloc[0]
        for cl in CLASSES:
            s.loc[g.index, f"{cl}_obs_exp"] = g[cl] / bg[cl] if bg[cl] > 0 else np.nan
    s.to_csv(f"{a.out}_features.tsv", sep="\t", index=False)
    hot = pd.concat(hots, ignore_index=True)
    hot.to_csv(f"{a.out}_hotspots.tsv", sep="\t", index=False)
    pd.set_option("display.width", 220)
    log(s[["gene_set", "feature", "n"] + list(CLASSES)].round(4).to_string(index=False))
    log(s[["gene_set", "feature"] + [f"{c}_obs_exp" for c in CLASSES]].round(2).to_string(index=False))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    first = a.gene_set[0].split("=", 1)[0]
    hh = hot[hot["gene_set"] == first]
    fig, axes = plt.subplots(len(HOT), 2, figsize=(12, 2.6 * len(HOT)), squeeze=False)
    for i, fam in enumerate(HOT):
        for j, ori in enumerate(("sense", "antisense")):
            ax = axes[i][j]
            for feat, col in (("acceptor", "tab:red"), ("donor", "tab:blue"), ("TSS", "tab:green"), ("TES", "tab:gray")):
                d = hh[(hh.rfamily == fam) & (hh.orientation == ori) & (hh.feature == feat)]
                if len(d):
                    ax.bar(d["cons_pos"], d["n"], width=3, color=col, alpha=0.7, label=feat)
            ax.set_title(f"{first}: {fam}, {ori} to the transcript", fontsize=9)
            ax.set_xlabel("consensus position")
            if j == 0:
                ax.set_ylabel("features")
            if ax.get_legend_handles_labels()[0]:
                ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(f"{a.out}_hotspots.png", dpi=130)
    log(f"wrote {a.out}_features.tsv, _hotspots.tsv, _hotspots.png")


if __name__ == "__main__":
    main()
