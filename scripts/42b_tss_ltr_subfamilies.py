#!/usr/bin/env python3
"""Patch 16. Positional clustering of lncRNA transcription starts in LTRs, within single subfamilies.

Input: <tss_ltr>_tss.tsv of step 42 (one row per start in an LTR record). For each named LTR group as published
(step 42 groups), for the group without internal-region records (names containing '-int'), and for every
subfamily with at least --min-starts starts in sense records, the share of sense starts in the most used 50-nt
consensus window is compared with the same statistic after placing each start at a random position of its own copy
(--n-perm permutations), exactly as in step 42.
Output: TSV with set, starts, sense starts, share sense, most used window, starts in it, share, expected share, P.
"""
import argparse

import numpy as np
import pandas as pd


def cons_pos(start, end, strand, cs, ce, pos):
    frac = (pos - start) / np.maximum(end - start - 1, 1)
    return np.round(np.where(strand == "+", cs + frac * (ce - cs), ce - frac * (ce - cs)))


def max_window(cpos, width):
    c = np.sort(cpos[np.isfinite(cpos)])
    if not len(c):
        return 0, np.nan
    j = np.searchsorted(c, c + width, side="left")
    k = np.argmax(j - np.arange(len(c)))
    mid = (c[k] + c[j[k] - 1]) / 2
    return int(j[k] - k), float(np.floor(mid - width / 2))


def test(s, width, n_perm, rng):
    obs, x0 = max_window(s["cons_pos"].to_numpy(float), width)
    st, en = s["copy_start"].to_numpy(), s["copy_end"].to_numpy()
    sd, cs, ce = s["copy_strand"].to_numpy(), s["cs"].to_numpy(), s["ce"].to_numpy()
    null = np.empty(n_perm)
    for k in range(n_perm):
        pos = st + (rng.random(len(s)) * (en - st)).astype(np.int64)
        null[k] = max_window(cons_pos(st, en, sd, cs, ce, pos), width)[0]
    return {"window": f"{x0:.0f}-{x0 + width:.0f}", "window_starts": obs, "window_share": obs / len(s),
            "window_share_expected": null.mean() / len(s), "window_p": (1 + (null >= obs).sum()) / (1 + n_perm)}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tss", required=True, help="step-42 <out>_tss.tsv")
    ap.add_argument("--window", type=int, default=50)
    ap.add_argument("--n-perm", type=int, default=1000)
    ap.add_argument("--min-starts", type=int, default=20)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)
    t = pd.read_csv(a.tss, sep="\t", low_memory=False)
    rows = []

    def add(label, level, d):
        s = d[d["orientation"] == "sense"]
        r = {"set": label, "level": level, "starts": len(d), "sense": len(s),
             "share_sense": len(s) / len(d) if len(d) else np.nan}
        if len(s) >= a.min_starts:
            r.update(test(s, a.window, a.n_perm, rng))
        rows.append(r)
    for g in [x for x in t["group"].unique() if x != "other LTR"]:
        d = t[t["group"] == g]
        add(g, "group as published", d)
        d2 = d[~d["subfamily"].astype(str).str.contains("-int", regex=False)]
        if len(d2) != len(d):
            add(g, "group without internal-region records", d2)
        for sf, dd in d2.groupby("subfamily"):
            if (dd["orientation"] == "sense").sum() >= a.min_starts:
                add(sf, f"subfamily of {g}", dd)
    out = pd.DataFrame(rows)
    out.to_csv(a.out, sep="\t", index=False, float_format="%.4g")
    print(out.to_string(index=False, float_format=lambda x: f"{x:.3g}"))


if __name__ == "__main__":
    main()
