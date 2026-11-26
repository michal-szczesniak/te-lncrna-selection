#!/usr/bin/env python3
"""Revision step 25b. Neutral baseline for the age of splice sites: are splice-site dinucleotides kept at the
orthologous position in other species more often than the same dinucleotide where it has no function?

Input: the per-site table of 25_splice_age.py (<out>_sites.tsv.gz; presence of each site in every species), the
step-3b mapping and copies (map_grouped.parquet, copies_grouped.tsv.gz), the positions of step 18/19 (to exclude
copies that host any splice site), hg38 and the same --species list as step 25.

Controls (one site = one GT or AG dinucleotide in transcript order):
  TE-derived sites   the same two consensus positions (grouped coordinates) in unused intronic or intergenic
                     copies of the same original subfamily, carrying the same dinucleotide in the same orientation
                     relative to the TE (so the control reads GT/AG on a 'virtual' transcript strand); up to
                     --n-ctrl per site. Created sites thus get controls with the same derived dinucleotide.
  non-TE and protein-coding sites   GT (donors) or AG (acceptors) on the transcript strand inside the same
                     intron, >= --flank-min nt from both ends, outside repeats, with the same CpG status as the
                     site when possible (a base is in a CpG if it or its neighbour forms CG); up to --n-ctrl.
For every control the orthologous dinucleotide is read exactly as in step 25 (present / absent / unaligned);
for TE controls the presence of the control copy itself is also tested.

Per set, subset and species: share of sites present (of aligned) and of their controls, the paired
difference site - control (sites with at least one aligned control; mean of the aligned controls of each site),
95% interval from --n-boot resamples of genes, and the ratio of absence rates (site / control). For TE-derived
sites also conditional on the copy being present in that species (site copy and control copies).
Positive difference = splice sites are kept more often than the same dinucleotide without function.

Outputs: <out>_baseline.tsv, <out>_controls.tsv.gz
"""
import argparse
import importlib.util
import os
import re

import numpy as np
import pandas as pd
import pysam

from ridl_utils import log


def load(name):
    here = os.path.dirname(os.path.abspath(__file__))
    spec = importlib.util.spec_from_file_location(name.replace(".py", ""), os.path.join(here, name))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


M25 = load("25_splice_age.py")
COMP = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}
RCMOTIF = {"GT": "AC", "AG": "CT"}  # genomic + strand string of a motif read on the - strand


def tx_dinuc(fa, chrom, p1, p2, strand):
    b1, b2 = fa.fetch(chrom, p1, p1 + 1).upper(), fa.fetch(chrom, p2, p2 + 1).upper()
    return (b1 + b2) if strand == "+" else (COMP.get(b1, "N") + COMP.get(b2, "N"))


def cpg_flag(fa, chrom, p1, p2):
    lo, hi = min(p1, p2), max(p1, p2)
    s = fa.fetch(chrom, max(lo - 1, 0), hi + 2).upper()
    return "CG" in s


def in_cpg(seq):
    """per base of seq: True if the base is part of a CpG (it and its neighbour form CG)."""
    a = np.frombuffer(seq.encode(), dtype="S1")
    cg = np.zeros(len(a), dtype=bool)
    if len(a) > 1:
        start = (a[:-1] == b"C") & (a[1:] == b"G")
        cg[:-1] |= start
        cg[1:] |= start
    return cg


def presence(lo, gen, chrom, p1, p2, strand, motif):
    m1, m2 = M25.lift_one(lo, chrom, p1), M25.lift_one(lo, chrom, p2)
    if m1 is None or m2 is None or m1[0] != m2[0] or m1[2] != m2[2] or abs(m1[1] - m2[1]) != 1:
        return "unaligned"
    flip = (m1[2] == "-") != (strand == "-")
    b1, b2 = gen.base(m1[0], m1[1]), gen.base(m2[0], m2[1])
    if flip:
        b1, b2 = COMP.get(b1, "N"), COMP.get(b2, "N")
    if b1 + b2 == motif:
        return "present"
    return "absent" if "N" not in b1 + b2 else "unaligned"


def copy_presence(lo, chrom, start, end, trim):
    L = int(end) - int(start) - 2 * trim
    if L < 20:
        return "short"
    e1 = M25.lift_one(lo, chrom, int(start) + trim)
    e2 = M25.lift_one(lo, chrom, int(end) - 1 - trim)
    ok = (e1 is not None and e2 is not None and e1[0] == e2[0] and e1[2] == e2[2]
          and 0.5 * L <= abs(e2[1] - e1[1]) <= 2 * L)
    return "present" if ok else "absent"


def te_controls(S, a, fa, rng):
    te = S[S["set"] == "lnc_TE_site"].copy()
    if te.empty:
        return pd.DataFrame()
    copies = pd.read_csv(a.copies, sep="\t", dtype={"chrom": str},
                         usecols=["copy_id", "chrom", "start", "end", "te_strand", "subfamily_orig", "compartment"])
    mp = pd.read_parquet(a.map, columns=["copy_id", "pos", "cons_pos"])
    mp = mp.merge(copies[["copy_id", "chrom"]], on="copy_id")
    # consensus positions (grouped coordinates) of both bases of every TE site
    mte = mp[mp["copy_id"].isin(set(te["copy_id"]))][["copy_id", "pos", "cons_pos"]]
    te = te.merge(mte.rename(columns={"pos": "pos1", "cons_pos": "c1"}), on=["copy_id", "pos1"], how="left")
    te = te.merge(mte.rename(columns={"pos": "pos2", "cons_pos": "c2"}), on=["copy_id", "pos2"], how="left")
    te = te.dropna(subset=["c1", "c2"])
    te["c1"], te["c2"] = te["c1"].astype(int), te["c2"].astype(int)
    # copies hosting any position of the splice test (sites or flanks) are not controls
    P = pd.read_csv(a.positions, sep="\t", usecols=["chrom", "pos", "set"], dtype={"chrom": str}, low_memory=False)
    P = P[~P["set"].isin(["lnc_TE_ctrl"])]
    pk = pd.MultiIndex.from_frame(P[["chrom", "pos"]].drop_duplicates())
    hosting = set(mp.loc[pd.MultiIndex.from_frame(mp[["chrom", "pos"]]).isin(pk), "copy_id"])
    pool_c = copies[copies["compartment"].isin(["intron", "intergenic"]) & ~copies["copy_id"].isin(hosting)
                    & copies["subfamily_orig"].isin(set(te["subfamily_orig"]))]
    need = pd.concat([te[["subfamily_orig", "c1"]].rename(columns={"c1": "cons_pos"}),
                      te[["subfamily_orig", "c2"]].rename(columns={"c2": "cons_pos"})]).drop_duplicates()
    pool = mp[mp["copy_id"].isin(set(pool_c["copy_id"]))].merge(
        pool_c[["copy_id", "te_strand", "subfamily_orig", "start", "end"]], on="copy_id").merge(need, on=["subfamily_orig", "cons_pos"])
    log(f"TE controls: {len(te)} sites, pool copies {len(pool_c)}, pool bases at needed positions {len(pool)}")
    by = {k: g for k, g in pool.groupby(["subfamily_orig", "cons_pos"])}
    rows = []
    for (sf, c1, c2), g in te.groupby(["subfamily_orig", "c1", "c2"]):
        a1, a2 = by.get((sf, c1)), by.get((sf, c2))
        if a1 is None or a2 is None:
            continue
        j = a1.merge(a2[["copy_id", "pos"]], on="copy_id", suffixes=("1", "2"))
        j = j[(j["pos1"] - j["pos2"]).abs() == 1]
        if j.empty:
            continue
        if len(j) > a.max_cand:
            j = j.iloc[rng.choice(len(j), a.max_cand, replace=False)]
        for r in g.itertuples(index=False):
            # virtual transcript strand of the control copy: same relation to the TE as the site
            vs = np.where(r.te_orientation == "sense", j["te_strand"], np.where(j["te_strand"] == "+", "-", "+"))
            jj = j.assign(vstrand=vs)
            order = ((jj["vstrand"] == "+") & (jj["pos2"] == jj["pos1"] + 1)) | ((jj["vstrand"] == "-") & (jj["pos2"] == jj["pos1"] - 1))
            jj = jj[order.to_numpy()]
            if jj.empty:
                continue
            ok = [tx_dinuc(fa, c, int(p1), int(p2), s) == r.motif
                  for c, p1, p2, s in zip(jj["chrom"], jj["pos1"], jj["pos2"], jj["vstrand"])]
            jj = jj[np.array(ok, dtype=bool)]
            if jj.empty:
                continue
            k = min(a.n_ctrl, len(jj))
            pick = jj.iloc[rng.choice(len(jj), k, replace=False)]
            for q in pick.itertuples(index=False):
                rows.append({"set": r.set, "site_id": r.site_id, "ctrl_kind": "unused copy", "chrom": q.chrom,
                             "pos1": int(q.pos1), "pos2": int(q.pos2), "tx_strand": q.vstrand, "motif": r.motif,
                             "copy_id": q.copy_id, "copy_start": int(q.start), "copy_end": int(q.end)})
    C = pd.DataFrame(rows)
    log(f"TE controls: {len(C)} for {C['site_id'].nunique() if len(C) else 0} of {len(te)} sites")
    return C


def intron_controls(S, a, fa, rng):
    m18 = load("18_splice_prepare.py")
    sub = S[S["set"].isin(["lnc_nonTE_site", "pc_site"])].copy()
    if sub.empty:
        return pd.DataFrame()
    rm = m18.read_rmsk(a.rmsk_out, set(sub["chrom"]))
    m = sub["site_id"].astype(str).str.extract(r"^(?P<ichrom>[^:]+):(?P<istart>\d+)-(?P<iend>\d+)(?P<istrand>[+-])[da]$")
    sub = pd.concat([sub, m], axis=1).dropna(subset=["istart"])
    sub["istart"], sub["iend"] = sub["istart"].astype(int), sub["iend"].astype(int)
    rows = []
    F = a.flank_min
    for r in sub.itertuples(index=False):
        lo, hi = r.istart + F, r.iend - F
        if hi - lo < 20:
            continue
        if hi - lo > a.max_window:
            lo = int(rng.integers(lo, hi - a.max_window))
            hi = lo + a.max_window
        seq = fa.fetch(r.chrom, lo, hi).upper()
        pat = r.motif if r.tx_strand == "+" else RCMOTIF[r.motif]
        idx = np.array([mm.start() for mm in re.finditer(f"(?={pat})", seq)], dtype=np.int64)
        if len(idx) == 0:
            continue
        cg = in_cpg(seq)
        g1 = lo + idx
        fl_all = cg[idx] | cg[idx + 1]
        if r.tx_strand == "+":
            p1, p2 = g1, g1 + 1
        else:
            p1, p2 = g1 + 1, g1
        rep = (m18.hit_rmsk(r.chrom, p1, rm) >= 0) | (m18.hit_rmsk(r.chrom, p2, rm) >= 0)
        p1, p2, fl = p1[~rep], p2[~rep], fl_all[~rep]
        if len(p1) == 0:
            continue
        want = cpg_flag(fa, r.chrom, int(r.pos1), int(r.pos2))
        same = np.flatnonzero(fl == want)
        cand = same if len(same) else np.arange(len(p1))
        pick = rng.choice(cand, min(a.n_ctrl, len(cand)), replace=False)
        for k in pick:
            rows.append({"set": r.set, "site_id": r.site_id, "ctrl_kind": "intron dinucleotide", "chrom": r.chrom,
                         "pos1": int(p1[k]), "pos2": int(p2[k]), "tx_strand": r.tx_strand, "motif": r.motif,
                         "copy_id": "", "copy_start": -1, "copy_end": -1, "cpg_matched": bool(fl[k] == want)})
    C = pd.DataFrame(rows)
    if len(C):
        log(f"intron controls: {len(C)} for {C['site_id'].nunique()} of {len(sub)} sites; "
            f"CpG status matched {C['cpg_matched'].mean():.1%}")
    return C


def boot_ci(d, w, cl, n_boot, rng):
    ids, inv = np.unique(cl, return_inverse=True)
    out = []
    for _ in range(n_boot):
        cnt = np.bincount(rng.integers(0, len(ids), len(ids)), minlength=len(ids))[inv].astype(float)
        ww = cnt * w
        out.append((d * ww).sum() / ww.sum() if ww.sum() > 0 else np.nan)
    out = np.array(out, dtype=float)
    return (np.nanpercentile(out, 2.5), np.nanpercentile(out, 97.5)) if np.isfinite(out).any() else (np.nan, np.nan)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sites", required=True, help="25_splice_age.py <out>_sites.tsv.gz")
    ap.add_argument("--positions", required=True, help="positions of step 18/19 (copies hosting any are excluded)")
    ap.add_argument("--map", required=True, help="step-3b map_grouped.parquet")
    ap.add_argument("--copies", required=True, help="step-3b copies_grouped.tsv.gz")
    ap.add_argument("--rmsk-out", required=True)
    ap.add_argument("--genome-fasta", required=True)
    ap.add_argument("--species", action="append", default=[], help="name=chain=genome, as in step 25")
    ap.add_argument("--n-ctrl", type=int, default=3)
    ap.add_argument("--max-cand", type=int, default=60, help="candidates examined per site (or consensus pair)")
    ap.add_argument("--flank-min", type=int, default=60)
    ap.add_argument("--max-window", type=int, default=20000)
    ap.add_argument("--copy-trim", type=int, default=20)
    ap.add_argument("--n-boot", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)
    S = pd.read_csv(a.sites, sep="\t", low_memory=False, dtype={"chrom": str})
    log(f"sites: {S.groupby('set').size().to_dict()}")
    fa = pysam.FastaFile(a.genome_fasta)
    C = pd.concat([te_controls(S, a, fa, rng), intron_controls(S, a, fa, rng)], ignore_index=True)
    if C.empty:
        raise SystemExit("no controls found")
    from pyliftover import LiftOver
    names = []
    for item in a.species:
        name, chain, gpath = item.split("=")
        names.append(name)
        log(f"{name}: loading {chain}")
        lo = LiftOver(chain)
        gen = M25.Genome(gpath)
        C[f"site_{name}"] = [presence(lo, gen, c, p1, p2, s, m) for c, p1, p2, s, m in
                             zip(C["chrom"], C["pos1"], C["pos2"], C["tx_strand"], C["motif"])]
        C[f"copy_{name}"] = [copy_presence(lo, c, s, e, a.copy_trim) if k == "unused copy" else ""
                             for c, s, e, k in zip(C["chrom"], C["copy_start"], C["copy_end"], C["ctrl_kind"])]
        log(f"  {name}: controls " + ", ".join(f"{k} {v}" for k, v in C[f"site_{name}"].value_counts().items()))
        del lo
    C.to_csv(f"{a.out}_controls.tsv.gz", sep="\t", index=False)

    rows = []
    for name in names:
        for cond in ([False, True]):
            for st in ("lnc_TE_site", "lnc_nonTE_site", "pc_site"):
                if cond and st != "lnc_TE_site":
                    continue
                s = S[S["set"] == st]
                c = C[C["set"] == st].copy()
                if s.empty or c.empty:
                    continue
                c_ok = c[f"site_{name}"] != "unaligned"
                if cond:
                    c_ok &= c[f"copy_{name}"] == "present"
                c = c[c_ok]
                cm = c.assign(pr=(c[f"site_{name}"] == "present").astype(float)).groupby("site_id").agg(
                    ctrl_present=("pr", "mean"), n_ctrl=("pr", "size"))
                s_ok = s[f"site_{name}"] != "unaligned"
                if cond:
                    s_ok &= s[f"copy_{name}"] == "present"
                ss = s[s_ok].set_index("site_id")
                ss = ss.join(cm, how="inner")
                ss["site_present"] = (ss[f"site_{name}"] == "present").astype(float)
                cols = ["site_type"] + (["rep_family", "ancestral", "te_orientation"] if st == "lnc_TE_site" else [])
                groups = [("all", ss)] + [(f"{col}={v}", g) for col in cols for v, g in ss.groupby(col)]
                for sub, g in groups:
                    if len(g) < 10:
                        continue
                    d = (g["site_present"] - g["ctrl_present"]).to_numpy()
                    lo_, hi_ = boot_ci(d, np.ones(len(d)), g["gene_id"].astype(str).to_numpy(), a.n_boot, rng)
                    ps, pc = g["site_present"].mean(), g["ctrl_present"].mean()
                    rows.append({"set": st, "species": name, "condition": "copy present" if cond else "all aligned",
                                 "subset": sub, "sites": len(g), "controls": int(g["n_ctrl"].sum()),
                                 "site_present": ps, "ctrl_present": pc, "diff": ps - pc, "diff_lo": lo_,
                                 "diff_hi": hi_, "absence_ratio": (1 - ps) / (1 - pc) if pc < 1 else np.nan})
    T = pd.DataFrame(rows)
    T.to_csv(f"{a.out}_baseline.tsv", sep="\t", index=False)
    pd.set_option("display.width", 250)
    log(T[T["subset"] == "all"].round(4).to_string(index=False))


if __name__ == "__main__":
    main()
