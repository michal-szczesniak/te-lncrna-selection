#!/usr/bin/env python3
"""Revision step 41 (patch 6). Restrict the exonic compartment of a step-4 annotation to the exons of the Ensembl
canonical transcript of each gene.

The main analysis defines lncRNA exons as the union of the exons of all annotated transcripts. Rarely used
alternative exons can dilute constraint that acts on the main isoform. This step keeps intronic and intergenic
bases unchanged and drops exonic bases that do not lie in an exon of the canonical transcript (GTF tag
"Ensembl_canonical") of the same gene. Exons of non-canonical isoforms are therefore removed from the analysis, not
moved to the intronic compartment (the intronic compartment already excludes every annotated exon). Units keep
their ids; an exonic unit that loses all of its bases disappears.

Inputs: a step-4 annotation directory (one parquet file per chromosome; columns chrom, pos (0-based), compartment,
gene_id, ...), the GTF used in step 1. Output: the same files in --out-dir, and <out-dir>/canonical_summary.tsv
(exonic bases and units per gene before and after, and genes without a canonical transcript).
"""
import argparse
import glob
import os

import numpy as np
import pandas as pd

from ridl_utils import log, merge_intervals, read_gtf


def canonical_exons(gtf, chroms):
    genes, tx, ex = read_gtf(gtf, chroms)
    can = set(tx.loc[tx["canonical"], "transcript_id"])
    e = ex[ex["transcript_id"].isin(can)][["chrom", "start", "end", "gene_id"]]
    e = merge_intervals(e, by=("gene_id", "chrom"))
    out = {}
    for (gid, chrom), d in e.groupby(["gene_id", "chrom"]):
        d = d.sort_values("start")
        out[gid] = (d["start"].to_numpy(np.int64), d["end"].to_numpy(np.int64))
    log(f"canonical transcripts: {len(can)}; genes with canonical exons: {len(out)}")
    return out


def inside(pos, iv):
    s, e = iv
    k = np.searchsorted(s, pos, side="right") - 1
    ok = k >= 0
    res = np.zeros(len(pos), dtype=bool)
    res[ok] = pos[ok] < e[k[ok]]
    return res


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--annot", required=True, help="step-4 annotation directory")
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--out-dir", required=True)
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    files = sorted(glob.glob(os.path.join(a.annot, "*.parquet")))
    if not files:
        raise SystemExit(f"no parquet files in {a.annot}")
    chroms = set()
    for f in files:
        c = pd.read_parquet(f, columns=["chrom"])["chrom"]
        if len(c):
            chroms.add(str(c.iloc[0]))
    can = canonical_exons(a.gtf, sorted(chroms))
    rows = []
    for f in files:
        d = pd.read_parquet(f)
        if d.empty:
            d.to_parquet(os.path.join(a.out_dir, os.path.basename(f)), index=False)
            continue
        keep = np.ones(len(d), dtype=bool)
        ex = (d["compartment"].astype(str) == "exon").to_numpy()
        gid = d["gene_id"].astype(str).to_numpy()
        pos = d["pos"].to_numpy(np.int64)
        idx_ex = np.flatnonzero(ex)
        for g, rows_g in pd.Series(idx_ex).groupby(gid[idx_ex]):
            r = rows_g.to_numpy()
            iv = can.get(g)
            keep[r] = inside(pos[r], iv) if iv is not None else False
        before = d[ex].groupby("gene_id", observed=True).agg(bases_before=("pos", "size"),
                                                             units_before=("unit_id", "nunique"))
        kept = d[keep]
        after = kept[kept["compartment"].astype(str) == "exon"].groupby("gene_id", observed=True).agg(
            bases_after=("pos", "size"), units_after=("unit_id", "nunique"))
        r = before.join(after, how="left").fillna(0).reset_index()
        r["has_canonical"] = r["gene_id"].astype(str).isin(can.keys())
        rows.append(r)
        kept.to_parquet(os.path.join(a.out_dir, os.path.basename(f)), index=False)
        log(f"  {os.path.basename(f)}: exonic bases {int(ex.sum())} -> {int((keep & ex).sum())}")
    if rows:
        s = pd.concat(rows, ignore_index=True)
        s.to_csv(os.path.join(a.out_dir, "canonical_summary.tsv"), sep="\t", index=False)
        tot_b, tot_a = s["bases_before"].sum(), s["bases_after"].sum()
        log(f"exonic bases kept: {tot_a:,.0f} of {tot_b:,.0f} ({100 * tot_a / max(tot_b, 1):.1f}%); genes with exonic "
            f"bases left: {(s['bases_after'] > 0).sum()} of {len(s)}; without canonical transcript: "
            f"{(~s['has_canonical']).sum()}")


if __name__ == "__main__":
    main()
