#!/usr/bin/env python3
"""Revision step 27. Splice-site constraint by junction usage and by site age (and with alternative controls).

Input: annotated positions from 19_splice_analyze.py (<out>_positions.tsv.gz), junction usage from
26_junction_usage.py and, optionally, site ages from 25_splice_age.py. The comparisons of step 19 are rerun
(same estimators and cluster bootstrap, imported from 19_splice_analyze.py) within strata:
  usage   quantile bins of --usage-metric computed within each site set (lncRNA non-TE, protein-coding,
          lncRNA TE-derived); unmatched / zero-read junctions form their own bin 'no reads'
  age     the age classes of step 25 (TE-derived and non-TE lncRNA sites, protein-coding sites)
  ancestral  preformed vs created TE-derived sites
Splice sites are compared with their own controls as in step 19: flank bases of the same introns (non-TE and
protein-coding sites) or the same consensus position and base in unused copies (TE-derived sites; only the
strata of the selected sites are used). If the usage of a site's flank control is needed, it is that of the
site's own intron.
With --positions-alt (positions annotated by step 19 after 18b_splice_controls.py) the TE comparison is also
run against those alternative controls (e.g. intronic copies of the same genes).
Robustness of the ancestral contrast: created and preformed TE-derived sites are also split by TE family, site
type and orientation (subsets with >= 20 sites), and all TE comparisons of 'all' and 'ancestral=...' are repeated
with the matching strata additionally split by CpG context of the base (comparison '... (strata x CpG)').

Output: <out>_strata.tsv
"""
import argparse
import importlib.util
import os

import numpy as np
import pandas as pd

from ridl_utils import log


def load_19():
    here = os.path.dirname(os.path.abspath(__file__))
    spec = importlib.util.spec_from_file_location("splice19", os.path.join(here, "19_splice_analyze.py"))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def intron_key(site_id):
    return site_id.str[:-1]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--positions", required=True)
    ap.add_argument("--usage", required=True, help="26_junction_usage.py output")
    ap.add_argument("--age", default=None, help="25_splice_age.py <out>_sites.tsv.gz")
    ap.add_argument("--positions-alt", action="append", default=[], help="label=path (step-19 positions after 18b)")
    ap.add_argument("--usage-metric", default="total_reads")
    ap.add_argument("--bins", type=int, default=3)
    ap.add_argument("--n-boot", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    m19 = load_19()
    rng = np.random.default_rng(a.seed)
    pos = pd.read_csv(a.positions, sep="\t", low_memory=False)
    metrics = [c for c in pos.columns if c.startswith("phylop_")]
    us = pd.read_csv(a.usage, sep="\t")
    umap = us.set_index("intron_key")[a.usage_metric]
    pos["_intron"] = intron_key(pos["site_id"].astype(str))
    # TE controls have stratum ids as site_id -> no usage; handled below
    pos["usage"] = pos["_intron"].map(umap)

    def bins_within(s):
        v = s.astype(float)
        out = pd.Series("no reads", index=s.index, dtype=object)
        ok = v > 0
        if ok.sum() >= a.bins * 10:
            q = pd.qcut(v[ok].rank(method="first"), a.bins, labels=[f"usage q{i + 1}" for i in range(a.bins)])
            out[ok] = q.astype(str)
        return out

    rows = []

    def add(name, subset, d, stratum, cluster):
        if d["lab"].nunique() < 2:
            return
        r = m19.compare(d, "lab", stratum, cluster, metrics, a.n_boot, rng)
        rows.append({"comparison": name, "subset": subset, **r})

    age = None
    if a.age:
        age = pd.read_csv(a.age, sep="\t", low_memory=False)
        age = age.set_index(["set", "site_id"])

    # ------------------------------------------------------------------ flank designs (non-TE lncRNA, protein-coding)
    for name, s_set, c_set in (("lncRNA non-TE splice sites vs intronic flank", "lnc_nonTE_site", "lnc_nonTE_flank"),
                               ("protein-coding splice sites vs intronic flank", "pc_site", "pc_flank")):
        d = pos[pos["set"].isin([s_set, c_set])].copy()
        if d.empty:
            continue
        d["lab"] = (d["set"] == s_set).astype(int)
        site_usage = d[d["lab"] == 1].drop_duplicates("site_id").set_index("site_id")["usage"]
        ub = bins_within(site_usage.fillna(0))
        d["ubin"] = d["site_id"].map(ub)
        for b, dd in d.groupby("ubin"):
            add(name, f"{b}", dd, "ctx", "site_id")
        if age is not None:
            ag = age.loc[s_set]["age_class"] if s_set in age.index.get_level_values(0) else None
            if ag is not None:
                d["age"] = d["site_id"].map(ag)
                for b, dd in d.groupby("age"):
                    add(name, f"age={b}", dd, "ctx", "site_id")

    # ------------------------------------------------------------------ TE-derived sites vs unused copies
    blocks = [("unused copies (main)", pos)]
    for item in a.positions_alt:
        lab, path = item.split("=", 1)
        alt = pd.read_csv(path, sep="\t", low_memory=False)
        alt["_intron"] = intron_key(alt["site_id"].astype(str))
        alt["usage"] = alt["_intron"].map(umap)
        blocks.append((lab, alt))
    at = None
    if age is not None and "lnc_TE_site" in age.index.get_level_values(0):
        at = age.loc["lnc_TE_site"]
    for clab, P in blocks:
        te = P[P["set"] == "lnc_TE_site"].copy()
        ctrl = P[P["set"] == "lnc_TE_ctrl"].copy()
        if te.empty or ctrl.empty:
            continue
        su = te.drop_duplicates("site_id").set_index("site_id")["usage"]
        te["ubin"] = te["site_id"].map(bins_within(su.fillna(0)))
        subsets = {"all": te.index}
        for b in sorted(te["ubin"].dropna().unique()):
            subsets[b] = te.index[te["ubin"] == b]
        if at is not None:
            te["age"] = te["site_id"].map(at["age_class"])
            te["ancestral"] = te["site_id"].map(at["ancestral"])
            for b in sorted(te["age"].dropna().unique()):
                subsets[f"age={b}"] = te.index[te["age"] == b]
            for b in sorted(te["ancestral"].dropna().unique()):
                subsets[f"ancestral={b}"] = te.index[te["ancestral"] == b]
            for b in sorted(te["ancestral"].dropna().unique()):
                for col in ("rep_family", "site_type", "te_orientation"):
                    if col not in te:
                        continue
                    for v in sorted(te[col].dropna().unique()):
                        idx = te.index[(te["ancestral"] == b) & (te[col] == v)]
                        if te.loc[idx, "site_id"].nunique() >= 20:
                            subsets[f"ancestral={b}&{col}={v}"] = idx
        for sub, idx in subsets.items():
            s = te.loc[idx]
            c = ctrl[ctrl["stratum"].isin(set(s["stratum"]))]
            if s.empty or c.empty:
                continue
            d = pd.concat([s.assign(lab=1), c.assign(lab=0)], ignore_index=True)
            add(f"lncRNA TE-derived splice sites vs {clab}", sub, d, "stratum", "stratum")
            if "cpg" in d and (sub == "all" or (sub.startswith("ancestral=") and "&" not in sub)):
                d["stratum_cpg"] = d["stratum"].astype(str) + "|" + d["cpg"].astype(str)
                add(f"lncRNA TE-derived splice sites vs {clab} (strata x CpG)", sub, d, "stratum_cpg", "stratum")
    st = pd.DataFrame(rows)
    st.to_csv(f"{a.out}_strata.tsv", sep="\t", index=False)
    pd.set_option("display.width", 250)
    show = ["comparison", "subset", "clusters", "bases_site", "poly_site", "poly_ctrl", "mh_log_or", "mh_log_or_lo",
            "mh_log_or_hi", "singleton_diff"] + [c for c in st.columns if c.endswith("_diff") and c != "singleton_diff"]
    log(st[[c for c in show if c in st.columns]].round(4).to_string(index=False))


if __name__ == "__main__":
    main()
