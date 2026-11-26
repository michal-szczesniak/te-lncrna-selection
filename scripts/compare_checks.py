"""Compare whole-element results across the main analysis and the robustness runs (missing runs are skipped)."""
import os

import numpy as np
import pandas as pd

pd.set_option("display.width", 230)
rng = np.random.default_rng(1)
GROUPS = ["AluS", "MIR", "L2abc", "AluJ"]
RUNS = [("main", "work/prof_v5/{g}_b20_paired", "work/comp/te_{g}_paired"),
        ("within gene", "work/prof_gene/{g}_b20", None),
        ("TSS 100 nt", "work/prof_tss100/{g}_b20_paired", "work/comp_tss100/te_{g}_paired"),
        ("extended", "work/prof_ext/{g}_b20_paired", "work/comp_ext/te_{g}_paired"),
        ("extended, within gene", "work/prof_ext/{g}_b20_gene", None)]
TRENDS = [("main", "work/trend/te_vs_nonte_genes.tsv", "work/comp/te_{g}_paired"),
          ("TSS 100 nt", "work/trend_tss100/te_vs_nonte_genes.tsv", "work/comp_tss100/te_{g}_paired"),
          ("extended", "work/trend_ext/te_vs_nonte_genes.tsv", "work/comp_ext/te_{g}_paired")]


def n_rows(p):
    try:
        return len(pd.read_csv(p, sep="\t"))
    except (FileNotFoundError, pd.errors.EmptyDataError):
        return 0


def gc_matched(prefix):
    f = f"{prefix}_summary.tsv"
    if not os.path.exists(f):
        return None
    return pd.read_csv(f, sep="\t").set_index(["metric", "comparison"]).loc[("phylop_p447", "gc_matched")]


def intron_rate(prefix):
    f = f"{prefix}_bins.tsv"
    if not os.path.exists(f):
        return np.nan
    b = pd.read_csv(f, sep="\t")
    return (b.poly_c_noncpg_intron * b.units_intron).sum() / b.units_intron.sum()


rows = []
for g in GROUPS:
    for tag, prof, comp in RUNS:
        f = os.path.join(prof.format(g=g), "whole_element.tsv")
        if os.path.exists(f):
            w = pd.read_csv(f, sep="\t")
            for _, r in w[w.layer.isin(["population", "phylop_p447"])].iterrows():
                rows.append({"group": g, "measure": f"{r.layer} {r.orientation}", "run": tag, "n": r.units_fg,
                             "estimate": r.effect, "ci_lo": r.ci_lo, "ci_hi": r.ci_hi})
            rows.append({"group": g, "measure": "centered peaks", "run": tag,
                         "estimate": n_rows(os.path.join(prof.format(g=g), "peaks_centered.tsv"))})
        s = gc_matched(comp.format(g=g)) if comp else None
        if s is not None:
            rows.append({"group": g, "measure": "GC-matched phyloP", "run": tag,
                         "estimate": s.exon_minus_intron, "ci_lo": s.ci_lo, "ci_hi": s.ci_hi})
for tag, f, comp in TRENDS:
    if not os.path.exists(f):
        continue
    gdf = pd.read_csv(f, sep="\t")
    for sf, s in gdf[gdf.orientation == "pooled"].groupby("subfamily"):
        g = sf.replace("*", "")
        x = s["excess_poly_noncpg"].dropna().to_numpy()
        b = np.array([rng.choice(x, len(x)).mean() for _ in range(2000)])
        lo, hi = np.percentile(b, [2.5, 97.5])
        rows.append({"group": g, "measure": "TE excess, variants", "run": tag, "n": len(x),
                     "estimate": x.mean(), "ci_lo": lo, "ci_hi": hi})
        rows.append({"group": g, "measure": "largest TE-specific drop (%)", "run": tag,
                     "estimate": -100 * lo / intron_rate(comp.format(g=g))})
for tag, prof, comp in (("main", "work/prof_nonte/exon_vs_intron_paired", "work/comp/lnc_paired"),
                        ("extended", "work/prof_nonte_ext/exon_vs_intron_paired", "work/comp_ext/lnc_paired")):
    f = os.path.join(prof, "whole_element.tsv")
    if os.path.exists(f):
        w = pd.read_csv(f, sep="\t")
        for _, r in w[w.layer.isin(["population", "phylop_p447"])].iterrows():
            rows.append({"group": "ordinary", "measure": f"{r.layer}", "run": tag, "n": r.units_fg,
                         "estimate": r.effect, "ci_lo": r.ci_lo, "ci_hi": r.ci_hi})
    s = gc_matched(comp)
    if s is not None:
        rows.append({"group": "ordinary", "measure": "GC-matched phyloP", "run": tag,
                     "estimate": s.exon_minus_intron, "ci_lo": s.ci_lo, "ci_hi": s.ci_hi})
t = pd.DataFrame(rows)
order = {r[0]: i for i, r in enumerate(RUNS)}
t["_o"] = t["run"].map(order)
t = t.sort_values(["group", "measure", "_o"]).drop(columns="_o")
print(t.round(4).to_string(index=False))
