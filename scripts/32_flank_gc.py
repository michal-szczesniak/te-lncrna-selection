#!/usr/bin/env python3
"""Revision step 32. GC matching on the flanks of each unit, as a robustness check of the GC-matched
conservation (Fig. 3a).

14_composition.py bins units (copies, copy parts or repeat-free windows) by the GC content of the unit itself.
Substitutions change that GC content and are also what phyloP measures, so matching on it partly conditions on
the outcome. Here every unit also gets the GC content of its genomic flanks (--flank nt on each side, the unit
itself excluded), which reflects the regional base composition but not the substitutions inside the unit.

Exon and intron units are compared
  unmatched     difference of means;
  gc_own        own-GC quintiles, intron means re-weighted to the exon distribution (= 14_composition.py);
  gc_flank      the same with flank-GC quintiles;
  gc_both       own-GC terciles x flank-GC terciles (9 strata).
Quantiles are taken over exon and intron units together. A stratum with exon but no intron units makes the
matched difference undefined (NaN), as in 14_composition.py. 95% intervals come from --n-boot resamples of genes.
Metrics: mean phyloP of each track (all bases and non-CpG bases), and the share of variable callable sites.

Outputs: <out>_summary.tsv (metric, comparison, exon_minus_intron, ci_lo, ci_hi),
         <out>_bins.tsv (flank-GC quintiles: range, units and means per compartment),
         <out>_units.tsv.gz (unit_id, gene_id, compartment, own and flank GC, metrics).
"""
import argparse
import glob
import os

import numpy as np
import pandas as pd
import pysam

from ridl_stats import load_annot
from ridl_utils import log


def flank_gc(units, fasta, flank, min_bases):
    """GC share of the flanks [lo - flank, lo) and [hi, hi + flank) of every unit (0-based, hi exclusive)."""
    fa = pysam.FastaFile(fasta)
    lens = dict(zip(fa.references, fa.lengths))
    out = np.full(len(units), np.nan)
    for i, (c, lo, hi) in enumerate(zip(units["chrom"], units["lo"], units["hi"])):
        if c not in lens:
            continue
        s = (fa.fetch(c, max(0, lo - flank), lo) + fa.fetch(c, hi, min(lens[c], hi + flank))).upper()
        gc = s.count("G") + s.count("C")
        at = s.count("A") + s.count("T")
        if gc + at >= min_bases:
            out[i] = gc / (gc + at)
    fa.close()
    return out


def qbins(x, k):
    """Equal-count bins 0..k-1 over the finite values of x (ties broken by order, as pd.qcut on ranks); -1 = NaN."""
    b = np.full(len(x), -1, dtype=int)
    ok = np.isfinite(x)
    if ok.sum() >= k:
        r = pd.Series(x[ok]).rank(method="first").to_numpy()
        b[ok] = pd.qcut(r, k, labels=False)
    return b


def matched_diff(v, is_exon, strata, n_strata, w):
    """sum over strata of exon share x (weighted exon mean - weighted intron mean); NaN if a stratum holding
    exon units has no intron units. v may contain NaN (ignored); w = bootstrap multiplicities."""
    ok = np.isfinite(v) & (strata >= 0) & (w > 0)
    s, e, vv, ww = strata[ok], is_exon[ok], v[ok], w[ok]
    key = s * 2 + e.astype(int)
    n = np.bincount(key, weights=ww, minlength=2 * n_strata).reshape(n_strata, 2)
    t = np.bincount(key, weights=ww * vv, minlength=2 * n_strata).reshape(n_strata, 2)
    nex = n[:, 1].sum()
    if nex == 0:
        return np.nan
    has_e = n[:, 1] > 0
    if np.any(has_e & (n[:, 0] == 0)):
        return np.nan
    me = t[has_e, 1] / n[has_e, 1]
    mi = t[has_e, 0] / n[has_e, 0]
    return float(np.sum(n[has_e, 1] / nex * (me - mi)))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--annot", required=True)
    ap.add_argument("--subfamilies", default="nonTE")
    ap.add_argument("--query", default=None)
    ap.add_argument("--paired-genes", action="store_true")
    ap.add_argument("--genome-fasta", required=True, help="indexed FASTA (.fai) of the reference genome")
    ap.add_argument("--flank", type=int, default=1000, help="flank length on each side of a unit (nt)")
    ap.add_argument("--min-flank-bases", type=int, default=200, help="minimum A/C/G/T bases in the flanks")
    ap.add_argument("--bins", type=int, default=5)
    ap.add_argument("--n-boot", type=int, default=200)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)

    import pyarrow.parquet as pq
    first = sorted(glob.glob(os.path.join(a.annot, "*.parquet")))[0] if os.path.isdir(a.annot) else a.annot
    schema = pq.read_schema(first).names
    pcols = [c for c in schema if c.startswith("phylop_")]
    cols = ["unit_id", "gene_id", "compartment", "subfamily", "chrom", "pos", "ref", "cpg"] + pcols + \
        [c for c in ("gn_poly", "gn_callable") if c in schema]
    df = load_annot(a.annot, columns=cols, query=a.query, subfamilies=a.subfamilies.split(","))
    df = df[df["compartment"].isin(["exon", "intron"])].copy()
    df["gene_id"] = df["gene_id"].astype(str)
    df["chrom"] = df["chrom"].astype(str)
    if a.paired_genes:
        g = df.groupby("gene_id")["compartment"].nunique()
        df = df[df["gene_id"].isin(g[g == 2].index)]
    if df.empty:
        raise SystemExit("no exon/intron units left -- check --subfamilies, --query and --paired-genes")
    df["is_gc"] = df["ref"].astype(str).str.upper().isin(["G", "C"])
    has_pop = "gn_poly" in df.columns
    if has_pop:
        df["poly_c"] = np.where(df["gn_callable"], df["gn_poly"].astype(float), np.nan)
        df["poly_c_noncpg"] = np.where(df["gn_callable"] & ~df["cpg"].astype(bool), df["gn_poly"].astype(float), np.nan)
    metrics = list(pcols) + (["poly_c", "poly_c_noncpg"] if has_pop else [])
    for c in pcols:
        df[f"{c}_noncpg"] = df[c].where(~df["cpg"].astype(bool))
    metrics += [f"{c}_noncpg" for c in pcols]

    keys = ["unit_id", "gene_id", "compartment"]
    u = df.groupby(keys, observed=True).agg(chrom=("chrom", "first"), lo=("pos", "min"), hi=("pos", "max"),
                                           gc=("is_gc", "mean"), **{m: (m, "mean") for m in metrics}).reset_index()
    u["hi"] = u["hi"] + 1
    del df
    log(f"units: exon {int((u.compartment == 'exon').sum())}, intron {int((u.compartment == 'intron').sum())}; "
        f"genes {u.gene_id.nunique()}; computing flank GC (+-{a.flank} nt)")
    u["flank_gc"] = flank_gc(u, a.genome_fasta, a.flank, a.min_flank_bases)
    ok = u["flank_gc"].notna()
    r = np.corrcoef(u.loc[ok, "gc"], u.loc[ok, "flank_gc"])[0, 1] if ok.sum() > 2 else np.nan
    for comp in ("exon", "intron"):
        s = u[u.compartment == comp]
        log(f"  {comp}: median own GC {s.gc.median():.3f}, median flank GC {s.flank_gc.median():.3f}")
    log(f"  correlation of own and flank GC: {r:.3f}; units without flank GC: {int((~ok).sum())}")

    is_exon = (u["compartment"] == "exon").to_numpy()
    b_own, b_flank = qbins(u["gc"].to_numpy(float), a.bins), qbins(u["flank_gc"].to_numpy(float), a.bins)
    t_own, t_flank = qbins(u["gc"].to_numpy(float), 3), qbins(u["flank_gc"].to_numpy(float), 3)
    both = np.where((t_own >= 0) & (t_flank >= 0), t_own * 3 + t_flank, -1)
    designs = {"unmatched": (np.zeros(len(u), dtype=int), 1), "gc_own": (b_own, a.bins),
               "gc_flank": (b_flank, a.bins), "gc_both": (both, 9)}
    # the unmatched and own-GC comparisons use every unit, as in 14_composition.py; flank designs drop units
    # without flank GC through their stratum -1
    vals = {m: u[m].to_numpy(float) for m in metrics}

    def stats(w):
        return {(m, d): matched_diff(vals[m], is_exon, s, k, w) for m in metrics for d, (s, k) in designs.items()}

    genes, gi = np.unique(u["gene_id"].to_numpy(), return_inverse=True)
    obs = stats(np.ones(len(u)))
    boots = []
    for _ in range(a.n_boot):
        cnt = np.bincount(rng.integers(0, len(genes), len(genes)), minlength=len(genes))
        boots.append(stats(cnt[gi].astype(float)))
    rows = []
    for (m, d), v in obs.items():
        bv = np.array([b[(m, d)] for b in boots], dtype=float)
        lo, hi = np.nanpercentile(bv, [2.5, 97.5]) if np.isfinite(bv).any() else (np.nan, np.nan)
        rows.append({"metric": m, "comparison": d, "exon_minus_intron": v, "ci_lo": lo, "ci_hi": hi,
                     "n_boot_finite": int(np.isfinite(bv).sum())})
    s = pd.DataFrame(rows)
    s.to_csv(f"{a.out}_summary.tsv", sep="\t", index=False)

    brows = []
    for b in range(a.bins):
        ub = u[b_flank == b]
        r_ = {"flank_gc_bin": b, "flank_gc_range": f"{ub['flank_gc'].min():.2f}-{ub['flank_gc'].max():.2f}",
              "units_exon": int((ub.compartment == "exon").sum()), "units_intron": int((ub.compartment == "intron").sum())}
        for m in metrics:
            r_[f"{m}_exon"] = ub.loc[ub.compartment == "exon", m].mean()
            r_[f"{m}_intron"] = ub.loc[ub.compartment == "intron", m].mean()
        brows.append(r_)
    pd.DataFrame(brows).to_csv(f"{a.out}_bins.tsv", sep="\t", index=False)
    u.drop(columns=["lo", "hi"]).to_csv(f"{a.out}_units.tsv.gz", sep="\t", index=False)
    pd.set_option("display.width", 200)
    show = s[s.metric.str.startswith("phylop") & ~s.metric.str.endswith("noncpg")]
    log(show.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
