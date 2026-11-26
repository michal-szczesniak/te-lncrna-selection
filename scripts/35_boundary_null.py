#!/usr/bin/env python3
"""Revision step 35 (patch 5). How often would an exonic TE copy contain an exon boundary by chance alone?

Fig. 4b reports the share of exonic copies whose span [start - pad, end + pad) contains a splice site of the host
gene (unit_class 'boundary' of step 28). A long copy that overlaps a short exon will often span a boundary just
because of its length. Null: the same copy (same length) placed uniformly at random so that it still overlaps the
same exon of the same gene (at least one base inside the exon core, i.e. the exon minus --core-pad nt at each end,
as the exonic compartment excludes the bases next to exon boundaries). The probability that the copy then contains
a splice site of the gene is computed exactly, by counting start positions, so the null has no sampling noise.
Orientation plays no role in the null, so the sense/antisense contrast of the observed shares is the part that
geometry cannot explain.

Input: the step-28 annotation subsets (exonic units with copy_id, gene_id, orientation, unit_class, pos), copy
coordinates (step 3b copies_grouped.tsv.gz), splice sites (step 1 splice_sites.tsv.gz: chrom, pos, gene_id) and the
GTF (exons of the host genes, merged per gene).

Output: <out>_units.tsv.gz (one row per exonic unit: observed class, null probability), <out>_summary.tsv (per
group x orientation: units, observed share 'boundary', expected share under the null, ratio, and 95% intervals of
both from resampling genes)
"""
import argparse
import glob
import os

import numpy as np
import pandas as pd

from ridl_utils import log, read_gtf


def merged_exons(ex):
    """gene_id -> sorted array of merged exon intervals (start, end), 0-based half-open."""
    out = {}
    for gid, g in ex.groupby("gene_id"):
        iv = g[["start", "end"]].to_numpy()
        iv = iv[np.argsort(iv[:, 0])]
        m = [list(iv[0])]
        for s, e in iv[1:]:
            if s <= m[-1][1]:
                m[-1][1] = max(m[-1][1], e)
            else:
                m.append([s, e])
        out[gid] = np.array(m, dtype=np.int64)
    return out


def null_prob(L, e0, e1, sites, pad, core_pad):
    """P(copy of length L placed at a uniform start x contains a splice site), x such that [x, x+L) overlaps the
    exon core [e0+core_pad, e1-core_pad) (the whole exon if the core is empty). A copy contains site p if
    x - pad <= p < x + L + pad, i.e. x in [p - L - pad + 1, p + pad]."""
    c0, c1 = e0 + core_pad, e1 - core_pad
    if c1 <= c0:
        c0, c1 = e0, e1
    lo, hi = c0 - L + 1, c1 - 1  # inclusive range of starts
    n = hi - lo + 1
    if n <= 0:
        return np.nan
    if len(sites) == 0:
        return 0.0
    a = sites - L - pad + 1
    b = sites + pad
    a, b = np.maximum(a, lo), np.minimum(b, hi)
    ok = b >= a
    if not ok.any():
        return 0.0
    a, b = a[ok], b[ok]
    o = np.argsort(a)
    a, b = a[o], b[o]
    covered, cur_a, cur_b = 0, a[0], b[0]
    for s, e in zip(a[1:], b[1:]):
        if s <= cur_b + 1:
            cur_b = max(cur_b, e)
        else:
            covered += cur_b - cur_a + 1
            cur_a, cur_b = s, e
    covered += cur_b - cur_a + 1
    return covered / n


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--annot-sub", action="append", required=True, help="GROUP=dir of a step-28 subset (repeatable)")
    ap.add_argument("--copies", required=True)
    ap.add_argument("--splice-sites", required=True)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--query", default=None)
    ap.add_argument("--pad", type=int, default=10, help="as --pad of step 28")
    ap.add_argument("--core-pad", type=int, default=10, help="bases next to exon ends excluded from the exonic compartment")
    ap.add_argument("--n-boot", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)
    cp = pd.read_csv(a.copies, sep="\t", usecols=["copy_id", "chrom", "start", "end"], dtype={"chrom": str})
    ss = pd.read_csv(a.splice_sites, sep="\t", dtype={"chrom": str})
    ss_by_gene = {g: np.sort(d["pos"].to_numpy().astype(np.int64)) for g, d in ss.groupby("gene_id")}

    units = []
    for item in a.annot_sub:
        grp, path = item.split("=", 1)
        files = sorted(glob.glob(os.path.join(path, "*.parquet"))) if os.path.isdir(path) else [path]
        for f in files:
            d = pd.read_parquet(f, columns=["unit_id", "copy_id", "compartment", "gene_id", "orientation",
                                            "unit_class", "chrom", "pos"], filters=[("compartment", "==", "exon")])
            d = d[d["compartment"] == "exon"]
            if a.query:
                d = d.query(a.query, engine="python")
            if d.empty:
                continue
            u = d.groupby(["unit_id", "copy_id", "gene_id", "orientation", "unit_class"], observed=True).agg(
                chrom=("chrom", "first"), lo=("pos", "min"), hi=("pos", "max")).reset_index()
            u["group"] = grp
            units.append(u)
    U = pd.concat(units, ignore_index=True)
    U = U.merge(cp.rename(columns={"chrom": "copy_chrom"}), on="copy_id", how="left")
    log(f"exonic units: {len(U)}; genes: {U.gene_id.nunique()}")
    chroms = sorted(U["chrom"].astype(str).unique())
    _, _, ex = read_gtf(a.gtf, chroms)
    ex = ex[ex["gene_id"].isin(set(U["gene_id"].astype(str)))]
    mex = merged_exons(ex)

    probs, obs_geo = [], []
    for r in U.itertuples(index=False):
        L = int(r.end) - int(r.start)
        exons = mex.get(str(r.gene_id))
        sites = ss_by_gene.get(str(r.gene_id), np.array([], dtype=np.int64))
        if exons is None or L <= 0:
            probs.append(np.nan)
            obs_geo.append(np.nan)
            continue
        # the exon that holds the unit's exonic bases
        k = np.flatnonzero((exons[:, 0] <= r.hi) & (exons[:, 1] > r.lo))
        if len(k) == 0:
            probs.append(np.nan)
            obs_geo.append(np.nan)
            continue
        e0, e1 = exons[k[0]]
        probs.append(null_prob(L, int(e0), int(e1), sites, a.pad, a.core_pad))
        s, e = int(r.start) - a.pad, int(r.end) + a.pad
        j = np.searchsorted(sites, s)
        obs_geo.append(float(j < len(sites) and sites[j] < e))
    U["p_null"] = probs
    U["boundary_obs"] = (U["unit_class"] == "boundary").astype(float)
    U["boundary_recomputed"] = obs_geo
    U = U[U["p_null"].notna()].copy()
    agree = (U["boundary_obs"] == U["boundary_recomputed"]).mean()
    log(f"units with a null: {len(U)}; observed class reproduced from the coordinates for {agree:.1%}")
    U.drop(columns=["lo", "hi"]).to_csv(f"{a.out}_units.tsv.gz", sep="\t", index=False)

    rows = []
    for (grp, ori), d in U.groupby(["group", "orientation"]):
        genes, gi = np.unique(d["gene_id"].astype(str).to_numpy(), return_inverse=True)
        o, n = d["boundary_obs"].to_numpy(), d["p_null"].to_numpy()
        bo, bn = [], []
        for _ in range(a.n_boot):
            w = np.bincount(rng.integers(0, len(genes), len(genes)), minlength=len(genes))[gi].astype(float)
            if w.sum() == 0:
                continue
            bo.append((o * w).sum() / w.sum())
            bn.append((n * w).sum() / w.sum())
        bo, bn = np.array(bo), np.array(bn)
        rows.append({"group": grp, "orientation": ori, "units": len(d), "genes": len(genes),
                     "observed": o.mean(), "obs_lo": np.percentile(bo, 2.5), "obs_hi": np.percentile(bo, 97.5),
                     "expected": n.mean(), "exp_lo": np.percentile(bn, 2.5), "exp_hi": np.percentile(bn, 97.5),
                     "ratio": o.mean() / n.mean() if n.mean() > 0 else np.nan,
                     "median_copy_len": float((d["end"] - d["start"]).median())})
    S = pd.DataFrame(rows)
    S.to_csv(f"{a.out}_summary.tsv", sep="\t", index=False)
    pd.set_option("display.width", 200)
    log(S.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
