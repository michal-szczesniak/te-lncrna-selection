#!/usr/bin/env python3
"""Revision step 33 (patch 5). Spike-in power scanned along the element, with genes resampled per replicate.

Same design and tests as step 22 (22_spikein_power.py, unchanged), with two additions:
  --scan LEN       plant segments of LEN nt starting at every --scan-step-th bin of the element (all testable
                   positions), instead of a few fixed segments: detection as a function of position.
  --boot-genes     every replicate resamples genes with replacement (one multiplicity per gene, shared by both
                   orientations); observed and permutation statistics of both orientations use these weights.
                   Without it (step 22) every replicate thins the same data, so detection is almost 0 or 1 and
                   r = 0 can only give 0 or 1. With it, replicates also vary the sampling of genes, and r = 0 gives
                   the share of resampled data sets in which the unmodified data already show a deficit there.
Thinning works as in step 22 (binomial removal of a fraction r of the variants in the planted bins of a random
fraction phi of the exonic units of one orientation).

Original description of step 22 follows.

Spike-in power on the real data: how strong and how long must a constrained segment be,
and in how many exonic copies, for the tests of the paper to detect it?

For one TE group, exonic (foreground) and intronic (background) units are taken exactly as in the main
same-gene profiles (06_profiles.py --paired-genes, strata tri x CpG x subfamily_orig, callable sites). A
constrained segment is planted by removing a fraction r of the gnomAD variants (binomial thinning of
polymorphic sites) in consensus bins [a, b] of a random fraction phi of the exonic units of one orientation.
The paper's tests are then rerun on the modified data:
  window     per-bin Mantel-Haenszel log OR, permutation p (unit labels), BH across all bins of both
             orientations; detected if a planted bin has q < --q and a negative effect
  centered   window effect minus the median window of the element (same permutations, BH); detected if a
             planted bin has q < --q and a negative centered effect
  whole      whole-element MH log OR, permutation p < --alpha and a negative effect
Repeating each setting --reps times gives detection rates (power). Thinning the aggregated
(unit x bin x stratum) counts is exactly equivalent to thinning variants base by base. Setting r = 0 gives
the false-positive rate of each test under the same design.

Segments are given in consensus positions and must start at a bin boundary (bin size --bin-size); 'all'
plants the reduction over the whole element (uniform constraint).

Outputs: <out>_power.tsv (one row per setting), <out>_reps.tsv (one row per replicate)
"""
import argparse

import numpy as np
import pandas as pd

from ridl_stats import Agg, load_annot, mh_logor, perm_pvalue, permuted_labels, select_comparison
from ridl_utils import bh_fdr, log


def strata_codes(d, cols):
    key = d[cols[0]].astype(str)
    for c in cols[1:]:
        key = key + "|" + d[c].astype(str)
    return pd.factorize(key)[0]


class Comparison:
    """One orientation: aggregated counts, observed and permuted per-bin and whole-element log ORs."""

    def __init__(self, d, nb, pop_strata, n_perm, scheme, min_units, rng):
        dp = d[d["gn_callable"]].copy()
        dp["_s"] = strata_codes(dp, pop_strata) if pop_strata else 0
        dp["gn_poly"] = dp["gn_poly"].astype(float)
        self.agg = Agg(dp, ["gn_poly"], strata_col="_s", n_bins=nb)
        self.base = self.agg.vals["gn_poly"].copy()
        self.nb = nb
        u1, u0 = self.agg.units_per_bin(self.agg.label)
        self.enough = (u1 >= min_units) & (u0 >= min_units)
        self.row_bin = self.agg.row_cell // self.agg.n_strata
        self.perms = [permuted_labels(self.agg, rng, scheme) for _ in range(n_perm)]
        self.ones = np.ones(len(self.agg.units))

    def stats(self, label, w=None):
        S = self.agg.sums(label, self.ones if w is None else w)
        a, c = S["gn_poly"]
        n1, n0 = S["_n"]
        return mh_logor(a, n1, c, n0, self.nb, self.agg.n_strata)

    def unit_weights(self, gene_mult):
        """Unit weights from gene multiplicities (dict gene -> count); units without a gene keep weight 1."""
        return np.array([gene_mult.get(c, 1.0) for c in self.agg.clusters], dtype=float)[self.agg.unit_cluster]

    def run(self, vals=None, w=None):
        """Observed and permutation statistics for the current values (vals replaces the variant counts);
        w = optional unit weights (gene bootstrap)."""
        self.agg.vals["gn_poly"] = self.base if vals is None else vals
        ob, oa = self.stats(self.agg.label, w)
        pb, pa = [], []
        for lab in self.perms:
            b, wa = self.stats(lab, w)
            pb.append(b)
            pa.append(wa)
        self.agg.vals["gn_poly"] = self.base
        pb, pa = np.array(pb), np.array(pa)
        with np.errstate(all="ignore"):  # as in step 6: minus the median of all windows of the element
            cen = ob - np.nanmedian(ob)
            pcen = pb - np.nanmedian(pb, axis=1)[:, None]
        res = dict(bin=np.where(self.enough, ob, np.nan), p=np.where(self.enough, perm_pvalue(ob, pb), np.nan),
                   cen=np.where(self.enough, cen, np.nan),
                   p_cen=np.where(self.enough, perm_pvalue(cen, pcen), np.nan),
                   all=oa, p_all=float(perm_pvalue(np.array([oa]), pa[:, None])[0]))
        return res

    def thinned(self, bins, r, phi, rng):
        """Variant counts after removing a fraction r of variants in `bins` of a fraction phi of fg units."""
        lab = self.agg.label
        chosen = (lab == 1) & (rng.random(len(lab)) < phi)
        rows = chosen[self.agg.row_unit] & np.isin(self.row_bin, bins)
        v = self.base.copy()
        if r > 0 and rows.any():
            v[rows] = rng.binomial(v[rows].astype(np.int64), 1.0 - r).astype(float)
        return v, int(chosen.sum())


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--annot", required=True)
    ap.add_argument("--subfamilies", required=True, help="one group, e.g. 'AluS*'")
    ap.add_argument("--orientation", default="sense", help="orientation that receives the planted segment")
    ap.add_argument("--segments", default="101-120,101-140,101-160,all",
                    help="comma list of consensus ranges a-b (bin-aligned) or 'all'")
    ap.add_argument("--reductions", default="0,0.1,0.25,0.5", help="fraction of variants removed")
    ap.add_argument("--fractions", default="1,0.5,0.25", help="fraction of exonic units carrying the segment")
    ap.add_argument("--reps", type=int, default=20)
    ap.add_argument("--n-perm", type=int, default=1000)
    ap.add_argument("--perm-scheme", choices=["global", "within_gene"], default="global")
    ap.add_argument("--bin-size", type=int, default=20)
    ap.add_argument("--min-units", type=int, default=20)
    ap.add_argument("--pop-strata", default="tri,cpg,subfamily_orig")
    ap.add_argument("--query", default=None)
    ap.add_argument("--no-paired-genes", action="store_true", help="use all genes (default: same genes only)")
    ap.add_argument("--q", type=float, default=0.1)
    ap.add_argument("--alpha", type=float, default=0.05)
    ap.add_argument("--scan", type=int, default=0, help="segment length (nt) planted at every position (0 = use --segments)")
    ap.add_argument("--scan-step", type=int, default=1, help="step between scanned segment starts (bins)")
    ap.add_argument("--boot-genes", action="store_true", help="resample genes in every replicate")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)
    strata = [c for c in a.pop_strata.split(",") if c]
    need = ["unit_id", "subfamily", "compartment", "orientation", "gene_id", "cons_pos", "gn_callable", "gn_poly"]
    df = load_annot(a.annot, columns=need + strata, query=a.query, subfamilies=[a.subfamilies])
    df = df[df["compartment"].isin(["exon", "intron"])].copy()
    df["bin"] = (df["cons_pos"] - 1) // a.bin_size
    nb = int(df["bin"].max()) + 1
    comps = {}
    for ori in ("sense", "antisense"):
        d = select_comparison(df, a.subfamilies, ori, "exon", "intron")
        if not a.no_paired_genes:
            g1 = set(d.loc[d["_lab"] == 1, "gene_id"]) - {"."}
            g0 = set(d.loc[d["_lab"] == 0, "gene_id"]) - {"."}
            d = d[d["gene_id"].isin(g1 & g0)]
        if d["_lab"].nunique() < 2:
            continue
        comps[ori] = Comparison(d, nb, strata, a.n_perm, a.perm_scheme, a.min_units, rng)
        log(f"{a.subfamilies} {ori}: fg units {(comps[ori].agg.label == 1).sum()}, "
            f"bg units {(comps[ori].agg.label == 0).sum()}, testable bins {int(comps[ori].enough.sum())}/{nb}")
    if a.orientation not in comps:
        raise SystemExit(f"no data for orientation {a.orientation}")
    other = [o for o in comps if o != a.orientation]
    base_other = {o: comps[o].run() for o in other}
    target = comps[a.orientation]
    genes_all = sorted(set().union(*[set(c.agg.clusters) for c in comps.values()]))

    def draw_weights():
        """Gene multiplicities shared by both orientations -> unit weights per comparison."""
        k = rng.integers(0, len(genes_all), len(genes_all))
        cnt = np.bincount(k, minlength=len(genes_all))
        mult = dict(zip(genes_all, cnt.astype(float)))
        return {o: comps[o].unit_weights(mult) for o in comps}

    segs = []
    seg_list = a.segments.split(",")
    if a.scan:
        if a.scan % a.bin_size:
            raise SystemExit(f"--scan {a.scan} is not a multiple of the bin size {a.bin_size}")
        k = a.scan // a.bin_size
        seg_list = [f"{b * a.bin_size + 1}-{(b + k) * a.bin_size}" for b in range(0, nb - k + 1, a.scan_step)]
    for s in seg_list:
        if s == "all":
            segs.append(("all", np.flatnonzero(target.enough)))
            continue
        lo, hi = (int(x) for x in s.split("-"))
        if (lo - 1) % a.bin_size or hi % a.bin_size:
            raise SystemExit(f"segment {s} is not aligned to {a.bin_size}-nt bins")
        bins = np.arange((lo - 1) // a.bin_size, hi // a.bin_size)
        bins = bins[(bins < nb)]
        bins = bins[target.enough[bins]]
        if len(bins) == 0:
            log(f"  segment {s}: no testable bins, skipped")
            continue
        segs.append((s, bins))
    reds = [float(x) for x in a.reductions.split(",")]
    fracs = [float(x) for x in a.fractions.split(",")]

    reps = []
    for seg, bins in segs:
        for r in reds:
            for phi in (fracs if r > 0 else [1.0]):
                for k in range(a.reps):
                    v, n_chosen = target.thinned(bins, r, phi, rng)
                    if a.boot_genes:
                        W = draw_weights()
                        res = target.run(v, W[a.orientation])
                        allres = {a.orientation: res, **{o: comps[o].run(None, W[o]) for o in other}}
                    else:
                        res = target.run(v)
                        allres = {a.orientation: res, **base_other}
                    # BH across all testable bins of both orientations, as in the main analysis
                    ps = np.concatenate([allres[o]["p"] for o in allres])
                    pcs = np.concatenate([allres[o]["p_cen"] for o in allres])
                    q = bh_fdr(ps)
                    qc = bh_fdr(pcs)
                    qt, qct = q[:nb], qc[:nb]  # the planted orientation comes first in allres
                    hit_w = bool(np.any((qt[bins] < a.q) & (res["bin"][bins] < 0)))
                    hit_c = bool(np.any((qct[bins] < a.q) & (res["cen"][bins] < 0)))
                    hit_a = bool((res["p_all"] < a.alpha) and (res["all"] < 0))
                    # false discoveries: bins outside the planted segment called in either direction
                    outside = np.setdiff1d(np.flatnonzero(target.enough), bins)
                    fd_w = int(np.sum(qt[outside] < a.q))
                    reps.append({"group": a.subfamilies, "orientation": a.orientation, "segment": seg,
                                 "seg_start": int(bins.min()) * a.bin_size + 1, "boot_genes": bool(a.boot_genes),
                                 "segment_nt": len(bins) * a.bin_size, "reduction": r, "fraction_units": phi,
                                 "rep": k, "units_planted": n_chosen,
                                 "effect_segment": float(np.nanmean(res["bin"][bins])),
                                 "effect_whole": float(res["all"]), "p_whole": res["p_all"],
                                 "min_q_segment": float(np.nanmin(qt[bins])),
                                 "detected_window": hit_w, "detected_centered": hit_c, "detected_whole": hit_a,
                                 "false_bins_window": fd_w})
                log(f"  {seg} r={r} phi={phi}: window {np.mean([x['detected_window'] for x in reps[-a.reps:]]):.2f}, "
                    f"centered {np.mean([x['detected_centered'] for x in reps[-a.reps:]]):.2f}, "
                    f"whole {np.mean([x['detected_whole'] for x in reps[-a.reps:]]):.2f}")
    R = pd.DataFrame(reps)
    R.to_csv(f"{a.out}_reps.tsv", sep="\t", index=False)
    keys = ["group", "orientation", "segment", "seg_start", "boot_genes", "segment_nt", "reduction", "fraction_units"]
    P = R.groupby(keys).agg(reps=("rep", "size"), effect_segment=("effect_segment", "mean"),
                            effect_whole=("effect_whole", "mean"), power_window=("detected_window", "mean"),
                            power_centered=("detected_centered", "mean"), power_whole=("detected_whole", "mean"),
                            false_bins_window=("false_bins_window", "mean")).reset_index()
    # expected reduction of the whole-element variant rate implied by the planted segment
    P["implied_whole_reduction"] = P["reduction"] * P["fraction_units"] * P["segment_nt"] / (nb * a.bin_size)
    P.to_csv(f"{a.out}_power.tsv", sep="\t", index=False)
    pd.set_option("display.width", 220)
    log(P.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
