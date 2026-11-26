#!/usr/bin/env python3
"""Revision step 42 (patch 6). Do LTR elements act as promoters of lncRNAs? Orientation and consensus position of
lncRNA transcription starts that lie in LTR retrotransposons.

An LTR that drives transcription of a lncRNA should (i) lie in the same orientation as the transcript and (ii) place
the start at a recurrent position of its consensus, where its own transcription start lies. A start that falls into
an LTR by chance should show neither. For every annotated transcript of the kept lncRNA genes (step-1
lnc_genes.tsv.gz, empty drop_reason) the first base is checked against RepeatMasker .out records of class LTR; the
subfamily, the orientation relative to the transcript and the consensus position (linear interpolation of the
RepeatMasker coordinates, as in step 17) are recorded. Starts are unique per gene and position.

Expectations:
  orientation   share of exonic sequence of the same genes that lies in sense LTRs among exonic sequence in LTRs
                (union of exons per gene, bases counted once)
  position      for every start in a sense LTR, a random position within the same LTR copy (uniform over the
                copy, the start's own copy); --n-perm permutations give the expected share of starts in the most
                used 50-nt consensus window of each subfamily group and an empirical P value

Groups (--groups, NAME:REGEX on the RepeatMasker subfamily name; the remaining LTR subfamilies form 'other LTR').

Outputs: <out>_tss.tsv (one row per start in an LTR), <out>_summary.tsv (per group: starts, share sense,
expected share sense, binomial P, most used window and its share of sense starts, expected share and P),
<out>_hist.tsv (counts per group, orientation and 10-nt consensus bin).
"""
import argparse
import re

import numpy as np
import pandas as pd
from scipy.stats import binomtest

from ridl_utils import log, merge_intervals, open_text, read_gtf

DEFAULT_GROUPS = "LTR12:^LTR12;THE1:^THE1;MLT1:^MLT1;MST:^MST;LTR7:^LTR7[A-Z]?$;MLT2:^MLT2"


def num(x):
    return int(x.strip("()"))


def read_ltr(path, chroms):
    rows = []
    with open_text(path) as fh:
        for line in fh:
            f = line.split()
            if len(f) < 15 or not f[0].isdigit() or f[4] not in chroms:
                continue
            cls = f[10].split("/")[0]
            if cls != "LTR":
                continue
            plus = f[8] == "+"
            if plus:
                cs, ce, clen = num(f[11]), num(f[12]), num(f[12]) + num(f[13])
            else:
                ce, cs, clen = num(f[12]), num(f[13]), num(f[12]) + num(f[11])
            rows.append((f[4], int(f[5]) - 1, int(f[6]), "+" if plus else "-", f[9], cs, ce, clen, int(f[14])))
    r = pd.DataFrame(rows, columns=["chrom", "start", "end", "strand", "subfamily", "cs", "ce", "clen", "rm_id"])
    return {c: g.sort_values("start").reset_index(drop=True) for c, g in r.groupby("chrom")}


def hit_records(chrom, pos, rm):
    """index of the record containing each position (-1: none)."""
    if chrom not in rm:
        return np.full(len(pos), -1)
    r = rm[chrom]
    S, E = r["start"].to_numpy(), r["end"].to_numpy()
    idx = np.searchsorted(S, pos, side="right") - 1
    hit = np.full(len(pos), -1)
    for k in range(8):
        cand = idx - k
        ok = (hit < 0) & (cand >= 0)
        ok[ok] &= E[cand[ok]] > pos[ok]
        hit[ok] = cand[ok]
    return hit


def cons_pos(rec, pos):
    frac = (pos - rec["start"].to_numpy()) / np.maximum(rec["end"].to_numpy() - rec["start"].to_numpy() - 1, 1)
    cs, ce = rec["cs"].to_numpy(), rec["ce"].to_numpy()
    return np.round(np.where(rec["strand"].to_numpy() == "+", cs + frac * (ce - cs), ce - frac * (ce - cs)))


def group_of(names, groups):
    out = np.full(len(names), "other LTR", dtype=object)
    for g, rx in reversed(groups):
        out[names.str.contains(rx, regex=True).to_numpy()] = g
    return out


def max_window(cpos, width):
    """largest number of positions in any window [x, x + width) and its start."""
    c = np.sort(cpos[np.isfinite(cpos)])
    if not len(c):
        return 0, np.nan
    j = np.searchsorted(c, c + width, side="left")
    k = np.argmax(j - np.arange(len(c)))
    # centre the window on the positions it contains
    mid = (c[k] + c[j[k] - 1]) / 2
    return int(j[k] - k), float(np.floor(mid - width / 2))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--lnc-genes", required=True, help="step-1 lnc_genes.tsv.gz")
    ap.add_argument("--rmsk-out", required=True)
    ap.add_argument("--groups", default=DEFAULT_GROUPS, help="';'-separated NAME:REGEX (first match wins)")
    ap.add_argument("--chroms", default=",".join(f"chr{i}" for i in range(1, 23)))
    ap.add_argument("--window", type=int, default=50)
    ap.add_argument("--n-perm", type=int, default=1000)
    ap.add_argument("--min-starts", type=int, default=20)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)
    chroms = [c for c in a.chroms.split(",") if c]
    groups = [tuple(x.split(":", 1)) for x in a.groups.split(";") if x]
    for _, rx in groups:
        re.compile(rx)

    lnc = pd.read_csv(a.lnc_genes, sep="\t", low_memory=False)
    kept = set(lnc.loc[lnc["drop_reason"].fillna("").astype(str) == "", "gene_id"])
    genes, tx, ex = read_gtf(a.gtf, set(chroms))
    tx = tx[tx["gene_id"].isin(kept)]
    ex = ex[ex["gene_id"].isin(kept)]
    plus = (tx["strand"] == "+").to_numpy()
    tss = pd.DataFrame({"gene_id": tx["gene_id"].to_numpy(), "chrom": tx["chrom"].to_numpy(),
                        "strand": tx["strand"].to_numpy(),
                        "pos": np.where(plus, tx["start"].to_numpy(), tx["end"].to_numpy() - 1)})
    tss = tss.drop_duplicates(["gene_id", "pos"]).reset_index(drop=True)
    log(f"kept genes {len(kept)}; transcripts {len(tx)}; unique starts {len(tss)}")
    rm = read_ltr(a.rmsk_out, set(chroms))

    # ------------------------------------------------------------- starts in LTRs
    rows = []
    for c, g in tss.groupby("chrom"):
        h = hit_records(c, g["pos"].to_numpy(), rm)
        m = h >= 0
        if not m.any():
            continue
        rec = rm[c].iloc[h[m]].reset_index(drop=True)
        gg = g[m].reset_index(drop=True)
        gg["subfamily"] = rec["subfamily"].to_numpy()
        gg["orientation"] = np.where(rec["strand"].to_numpy() == gg["strand"].to_numpy(), "sense", "antisense")
        gg["cons_pos"] = cons_pos(rec, gg["pos"].to_numpy())
        gg["cons_len"] = rec["clen"].to_numpy()
        gg["copy_start"], gg["copy_end"] = rec["start"].to_numpy(), rec["end"].to_numpy()
        gg["copy_strand"], gg["cs"], gg["ce"] = rec["strand"].to_numpy(), rec["cs"].to_numpy(), rec["ce"].to_numpy()
        rows.append(gg)
    t = pd.concat(rows, ignore_index=True)
    t["group"] = group_of(t["subfamily"].astype(str), groups)
    t.to_csv(f"{a.out}_tss.tsv", sep="\t", index=False)
    log(f"starts in LTRs: {len(t)} of {len(tss)} ({100 * len(t) / len(tss):.1f}%)")

    # ------------------------------------------------------------- expected orientation from exonic sequence
    ex_u = merge_intervals(ex[["chrom", "start", "end", "gene_id"]], by=("gene_id", "chrom"))
    strand = tx.drop_duplicates("gene_id").set_index("gene_id")["strand"]
    ex_u["strand"] = ex_u["gene_id"].map(strand)
    exp = {}
    for c, g in ex_u.groupby("chrom"):
        if c not in rm:
            continue
        r = rm[c]
        S, E = r["start"].to_numpy(), r["end"].to_numpy()
        grp = group_of(r["subfamily"].astype(str), groups)
        for s, e, st in zip(g["start"].to_numpy(), g["end"].to_numpy(), g["strand"].to_numpy()):
            i0 = max(np.searchsorted(S, s, side="right") - 8, 0)
            i1 = np.searchsorted(S, e, side="left")
            for i in range(i0, i1):
                ov = min(e, E[i]) - max(s, S[i])
                if ov > 0:
                    key = (grp[i], "sense" if r["strand"].iat[i] == st else "antisense")
                    exp[key] = exp.get(key, 0) + ov
    for key in list(exp):
        exp[("all LTR", key[1])] = exp.get(("all LTR", key[1]), 0) + exp[key]

    # ------------------------------------------------------------- summary per group
    out, hist = [], []
    t2 = pd.concat([t, t.assign(group="all LTR")], ignore_index=True)
    for grp, d in t2.groupby("group"):
        n, ns = len(d), int((d["orientation"] == "sense").sum())
        es, ea = exp.get((grp, "sense"), 0), exp.get((grp, "antisense"), 0)
        p0 = es / (es + ea) if es + ea else np.nan
        pb = binomtest(ns, n, p0).pvalue if n and np.isfinite(p0) else np.nan
        row = {"group": grp, "starts": n, "sense": ns, "share_sense": ns / n if n else np.nan,
               "exonic_bp_sense": es, "exonic_bp_antisense": ea, "expected_share_sense": p0, "p_binom": pb,
               "subfamilies": d["subfamily"].nunique()}
        s = d[d["orientation"] == "sense"]
        if grp != "all LTR" and len(s) >= a.min_starts:
            # the most used window on the consensus of the dominant subfamily length scale
            obs, x0 = max_window(s["cons_pos"].to_numpy(float), a.window)
            null = np.empty(a.n_perm)
            L = (s["copy_end"] - s["copy_start"]).to_numpy()
            for k in range(a.n_perm):
                pos = s["copy_start"].to_numpy() + (rng.random(len(s)) * L).astype(np.int64)
                rec = s[["copy_start", "copy_end", "copy_strand", "cs", "ce"]].rename(
                    columns={"copy_start": "start", "copy_end": "end", "copy_strand": "strand"})
                null[k] = max_window(cons_pos(rec, pos), a.window)[0]
            row.update({"window_start": x0, "window_end": x0 + a.window, "window_starts": obs,
                        "window_share": obs / len(s), "window_share_expected": null.mean() / len(s),
                        "window_p": (1 + (null >= obs).sum()) / (1 + a.n_perm),
                        "top_subfamily": s["subfamily"].value_counts().index[0]})
        out.append(row)
        for ori, dd in d.groupby("orientation"):
            b = (dd["cons_pos"] // 10 * 10).dropna().astype(int).value_counts().sort_index()
            hist.append(pd.DataFrame({"group": grp, "orientation": ori, "cons_bin": b.index, "starts": b.to_numpy()}))
    S = pd.DataFrame(out).sort_values("starts", ascending=False)
    S.to_csv(f"{a.out}_summary.tsv", sep="\t", index=False, float_format="%.4g")
    pd.concat(hist, ignore_index=True).to_csv(f"{a.out}_hist.tsv", sep="\t", index=False)
    log("\n" + S.to_string(index=False))


if __name__ == "__main__":
    main()
