#!/usr/bin/env python3
"""Sanity check against the literature: phyloP of lncRNA exons and introns, INCLUDING the regions our main
design excludes (TSS +-500 nt, splice-site margins), split by region and distance from the TSS.

For a random sample of genes, the reference transcript (step-1 ref_tx table) is split into: first exon,
internal exons, last exon (single-exon genes: 'single exon'), first intron, other introns. Repeat-masked
bases (all RepeatMasker records) are dropped, so this is ordinary sequence only. Exonic bases get their
distance from the TSS along the spliced transcript, intronic bases along the genome. Intronic sequence is
sampled (--intron-windows windows of --intron-window-bp per intron) to keep the run short.

Output: TSV with region, distance bin, bases, mean of each phyloP track and share of bases with phyloP > 2.
"""
import argparse

import numpy as np
import pandas as pd
import pyBigWig

from ridl_utils import index_by_chrom, log, merge_intervals, open_text

BINS = [0, 250, 500, 1000, 2000, 5000, 10 ** 9]
LABELS = ["0-250", "250-500", "500-1k", "1k-2k", "2k-5k", ">5k"]


def unmasked(chrom, s, e, ridx):
    """Boolean mask of positions s..e-1 NOT covered by repeats."""
    keep = np.ones(e - s, dtype=bool)
    if chrom in ridx:
        rs, re_ = ridx[chrom]
        i = np.searchsorted(re_, s, side="right")
        while i < len(rs) and rs[i] < e:
            keep[max(rs[i], s) - s:min(re_[i], e) - s] = False
            i += 1
    return keep


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ref-tx", required=True)
    ap.add_argument("--rmsk-out", required=True)
    ap.add_argument("--phylop", action="append", required=True, help="name=path.bw")
    ap.add_argument("--max-genes", type=int, default=4000)
    ap.add_argument("--intron-windows", type=int, default=20)
    ap.add_argument("--intron-window-bp", type=int, default=100)
    ap.add_argument("--chroms", default=",".join(f"chr{i}" for i in range(1, 23)))
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)
    chroms = set(a.chroms.split(","))

    tx = pd.read_csv(a.ref_tx, sep="\t", dtype={"chrom": str})
    tx = tx[tx["chrom"].isin(chroms)]
    genes = tx["gene_id"].unique()
    if len(genes) > a.max_genes:
        genes = rng.choice(genes, a.max_genes, replace=False)
    tx = tx[tx["gene_id"].isin(set(genes))]
    log(f"genes sampled: {len(genes)}")

    rows = []
    with open_text(a.rmsk_out) as fh:
        for line in fh:
            f = line.split()
            if len(f) >= 15 and f[0].isdigit() and f[4] in chroms:
                rows.append((f[4], int(f[5]) - 1, int(f[6])))
    ridx = index_by_chrom(merge_intervals(pd.DataFrame(rows, columns=["chrom", "start", "end"])))
    bws = {}
    for item in a.phylop:
        n, p = item.split("=", 1)
        bws[n] = pyBigWig.open(p)

    recs = []
    for gid, g in tx.groupby("gene_id"):
        g = g.sort_values("tx_offset")
        chrom, strand = g["chrom"].iloc[0], g["strand"].iloc[0]
        ex = g[["start", "end", "tx_offset"]].to_numpy()
        n_ex = len(ex)
        tss = ex[0][0] if strand == "+" else ex[0][1] - 1
        for k, (s, e, off) in enumerate(ex):
            region = "single exon" if n_ex == 1 else ("first exon" if k == 0 else
                                                      ("last exon" if k == n_ex - 1 else "internal exon"))
            pos = np.arange(s, e)
            dist = off + (pos - s if strand == "+" else e - 1 - pos)
            recs.append((chrom, s, e, region, dist))
        # introns between consecutive exons in transcript order
        for k in range(n_ex - 1):
            if strand == "+":
                s, e = ex[k][1], ex[k + 1][0]
            else:
                s, e = ex[k + 1][1], ex[k][0]
            if e - s < a.intron_window_bp:
                continue
            region = "first intron" if k == 0 else "other intron"
            W = a.intron_window_bp
            starts = rng.integers(s, e - W + 1, size=min(a.intron_windows, max(1, (e - s) // W)))
            for ws in np.unique(starts):
                pos = np.arange(ws, ws + W)
                recs.append((chrom, ws, ws + W, region, np.abs(pos - tss)))

    out = []
    for chrom, s, e, region, dist in recs:
        keep = unmasked(chrom, s, e, ridx)
        if not keep.any():
            continue
        vals = {}
        for n, bw in bws.items():
            try:
                vals[n] = np.array(bw.values(chrom, int(s), int(e), numpy=True), dtype=float)[keep]
            except RuntimeError:
                vals[n] = np.full(keep.sum(), np.nan)
        d = pd.DataFrame({"region": region, "dist": dist[keep], **vals})
        out.append(d)
    df = pd.concat(out, ignore_index=True)
    df["dist_bin"] = pd.cut(df["dist"], BINS, labels=LABELS, right=False)
    res = []
    for (region, b), g in df.groupby(["region", "dist_bin"], observed=True):
        r = {"region": region, "dist_from_tss": b, "bases": len(g)}
        for n in bws:
            v = g[n].dropna()
            r[f"mean_{n}"] = v.mean()
            r[f"share_{n}_gt2"] = (v > 2).mean()
        res.append(r)
    for region, g in df.groupby("region"):
        r = {"region": region, "dist_from_tss": "all", "bases": len(g)}
        for n in bws:
            v = g[n].dropna()
            r[f"mean_{n}"] = v.mean()
            r[f"share_{n}_gt2"] = (v > 2).mean()
        res.append(r)
    res = pd.DataFrame(res)
    res.to_csv(a.out, sep="\t", index=False)
    pd.set_option("display.width", 200)
    log(res.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
