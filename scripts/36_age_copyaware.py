#!/usr/bin/env python3
"""Revision step 36 (patch 5). Age of TE-derived splice sites when presence requires the copy itself.

Step 25 calls a site 'present' in another species when the lifted dinucleotide is GT/AG. For TE-derived sites this
can happen in a species that lacks the insertion, when the chain aligns non-orthologous sequence or a related
element (e.g. Alu-derived sites 'present' in mouse). Here a TE-derived site counts as present in a species only if
the copy is also present there (step 25 'copy_<species>' = present: both ends of the copy lift to the same
chromosome and strand at a plausible distance). A site that is GT/AG but whose copy is absent is set to
'no shared copy'. With --rule absent (default) this applies only when the copy was tested and found absent;
copies too short to test after trimming ('short') keep the site call. With --rule strict every copy not called
present (absent or short) is set to 'no shared copy'. Age classes are then recomputed with the rule of step 25
(deepest species with a present site; 'human-specific' if absent but aligned in the two closest species;
'unresolved' otherwise). Non-TE sets are unchanged. The summary also counts, per species and TE family, the
sites reclassified.

Input: 25_splice_age.py <out>_sites.tsv.gz. Outputs: <out>_sites.tsv.gz (added columns site_ca_<species>,
age_class_ca), <out>_summary.tsv (same layout as step 25: set, subset, sites, measure, n, share, lo, hi),
<out>_reclassified.tsv (per species and family: sites reclassified, of them with copy 'absent' and 'short')
"""
import argparse

import numpy as np
import pandas as pd

from ridl_utils import log


def wilson(k, n, z=1.96):
    if n == 0:
        return np.nan, np.nan, np.nan
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return p, c - h, c + h


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sites", required=True)
    ap.add_argument("--species", required=True, help="comma list, closest first (as in step 25)")
    ap.add_argument("--rule", default="absent", choices=["absent", "strict"],
                    help="absent: only copies tested and absent veto a site; strict: also copies too short to test")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    names = a.species.split(",")
    S = pd.read_csv(a.sites, sep="\t", low_memory=False)
    te = (S["set"] == "lnc_TE_site").to_numpy()
    rec = []
    for n in names:
        st = S[f"site_{n}"].astype(str).to_numpy()
        cp = S[f"copy_{n}"].astype(str).to_numpy() if f"copy_{n}" in S else np.full(len(S), "present")
        new = st.copy()
        veto = (cp == "absent") if a.rule == "absent" else (cp != "present")
        bad = te & (st == "present") & veto
        new[bad] = "no shared copy"
        S[f"site_ca_{n}"] = new
        for fam, k in pd.Series(S.loc[bad, "rep_family"]).value_counts().items():
            fb = bad & (S["rep_family"] == fam).to_numpy()
            rec.append({"species": n, "rep_family": fam, "reclassified": int(k),
                        "copy_absent": int((fb & (cp == "absent")).sum()), "copy_short": int((fb & (cp == "short")).sum()),
                        "te_sites_present_before": int((te & (st == "present") & (S["rep_family"] == fam)).sum())})
    pres = np.stack([S[f"site_ca_{n}"].to_numpy() == "present" for n in names], axis=1)
    alig = np.stack([~np.isin(S[f"site_ca_{n}"].to_numpy(), ["unaligned", "no shared copy"]) for n in names], axis=1)
    deepest = np.where(pres.any(axis=1), pres.shape[1] - 1 - np.argmax(pres[:, ::-1], axis=1), -1)
    age = np.where(deepest >= 0, np.array([f"shared with {n}" for n in names], dtype=object)[np.maximum(deepest, 0)],
                   "unresolved")
    if len(names) >= 2:
        age = np.where((deepest < 0) & alig[:, 0] & alig[:, 1], "human-specific", age)
    S["age_class_ca"] = np.where(te, age, S["age_class"].astype(str))
    S.to_csv(f"{a.out}_sites.tsv.gz", sep="\t", index=False)

    rows = []
    for stn, s in S.groupby("set"):
        subsets = [("all", s)]
        if stn == "lnc_TE_site":
            for col in ("rep_family", "te_orientation", "ancestral", "site_type"):
                if col in s:
                    subsets += [(f"{col}={v}", x) for v, x in s.groupby(col)]
        for sub, x in subsets:
            base = {"set": stn, "subset": sub, "sites": len(x)}
            for n in names:
                col = f"site_ca_{n}"
                al = ~x[col].isin(["unaligned", "no shared copy"])
                k, m = int((x[col] == "present").sum()), int(al.sum())
                p, lo, hi = wilson(k, m)
                rows.append({**base, "measure": f"present_in_{n} (of aligned, copy shared)", "n": m, "share": p,
                             "lo": lo, "hi": hi})
            for cl, k in x["age_class_ca"].value_counts().items():
                p, lo, hi = wilson(int(k), len(x))
                rows.append({**base, "measure": f"age={cl}", "n": len(x), "share": p, "lo": lo, "hi": hi})
    T = pd.DataFrame(rows)
    T.to_csv(f"{a.out}_summary.tsv", sep="\t", index=False)
    R = pd.DataFrame(rec)
    R.to_csv(f"{a.out}_reclassified.tsv", sep="\t", index=False)
    pd.set_option("display.width", 200)
    log(R.to_string(index=False) if len(R) else "no sites reclassified")
    log(T[T["subset"].isin(["all", "rep_family=Alu"]) & T["measure"].str.startswith("age=")].round(4).to_string(index=False))


if __name__ == "__main__":
    main()
