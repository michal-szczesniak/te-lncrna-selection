#!/usr/bin/env python3
"""Revision step 25. How old are lncRNA splice sites, and were the TE-derived ones already present in the TE?

Input: annotated splice-site positions from 19_splice_analyze.py (<out>_positions.tsv.gz): sets lnc_TE_site
(GT/AG dinucleotide inside a mapped TE copy), lnc_nonTE_site (outside any repeat) and pc_site (protein-coding).
One row per site is built from its two dinucleotide bases, in transcript order.

(1) Ancestral state of TE-derived sites. The consensus base of the copy's ORIGINAL subfamily (step 3
    map_align.parquet + consensus_rec.fa) is read at both bases and put in transcript orientation:
      preformed  the consensus already has GT (donor) or AG (acceptor) there: the TE brought the site
      created    at least one base differs from the consensus: the site arose by mutation after insertion
                 (or the consensus differs from the source element at that base)
(2) Presence in other species. Both bases are lifted with UCSC chain files (--species name=chain=genome;
    genome = .2bit or indexed FASTA) and read in transcript orientation: 'present' if the orthologous
    dinucleotide is GT/AG, 'absent' if both bases align but the dinucleotide differs, 'unaligned' otherwise.
    For TE-derived sites the copy itself is also lifted (two points --copy-trim nt inside its ends, same
    chromosome and strand, distance within 0.5-2x the copy length): 'copy present' means the insertion is
    shared with that species.
(3) Age class. Species must be listed from the closest to the most distant. The site is at least as old as the
    split from the most distant species in which it is present; 'human-specific' if absent but aligned in the
    two closest species; 'unresolved' if it is present nowhere and not aligned in both.
(4) Human polymorphism at the site: highest allele frequency (gnomAD AC/AN, callable bases) of a variant at
    either base; 'segregating' if any variant, 'common' if AF >= 1%.

Outputs: <out>_sites.tsv.gz (one row per site), <out>_summary.tsv (shares with Wilson 95% intervals)
"""
import argparse

import numpy as np
import pandas as pd
import pysam

from ridl_utils import log

COMP = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}
MOTIF = {"donor": "GT", "acceptor": "AG"}


def comp(b):
    return COMP.get(b, "N")


class Genome:
    def __init__(self, path):
        self.path = path
        if path.endswith(".2bit"):
            import py2bit
            self.tb = py2bit.open(path)
            self.fa = None
            self.chroms = set(self.tb.chroms().keys())
        else:
            self.fa = pysam.FastaFile(path)
            self.tb = None
            self.chroms = set(self.fa.references)

    def base(self, chrom, pos):
        if chrom not in self.chroms:
            return "N"
        if self.tb is not None:
            return self.tb.sequence(chrom, int(pos), int(pos) + 1).upper()
        return self.fa.fetch(chrom, int(pos), int(pos) + 1).upper()


def load_consensus(path):
    seqs = {}
    with pysam.FastxFile(path) as fh:
        for rec in fh:
            seqs[rec.name.split("#")[0]] = rec.sequence.upper()
    return seqs


def wilson(k, n, z=1.96):
    if n == 0:
        return np.nan, np.nan, np.nan
    p = k / n
    den = 1 + z * z / n
    c = (p + z * z / (2 * n)) / den
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return p, c - h, c + h


def lift_one(lo, chrom, pos):
    r = lo.convert_coordinate(chrom, int(pos))
    if not r:
        return None
    return r[0][0], int(r[0][1]), r[0][2]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--positions", required=True, help="19_splice_analyze.py <out>_positions.tsv.gz")
    ap.add_argument("--map-align", required=True, help="step-3 map_align.parquet (original consensus positions)")
    ap.add_argument("--consensus", required=True, help="step-3 consensus_rec.fa")
    ap.add_argument("--copies", required=True, help="step-3b copies_grouped.tsv.gz")
    ap.add_argument("--genome-fasta", required=True, help="hg38 FASTA (sanity check of the dinucleotides)")
    ap.add_argument("--species", action="append", default=[],
                    help="name=chain=genome, closest first (e.g. panTro6=hg38ToPanTro6.over.chain.gz=panTro6.2bit)")
    ap.add_argument("--copy-trim", type=int, default=20)
    ap.add_argument("--sets", default="lnc_TE_site,lnc_nonTE_site,pc_site")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    pos = pd.read_csv(a.positions, sep="\t", low_memory=False)
    pos = pos[pos["set"].isin(a.sets.split(","))].copy()
    pos = pos.sort_values(["set", "site_id", "dinuc_idx"])
    log(f"positions: {len(pos)} bases, {pos.groupby('set')['site_id'].nunique().to_dict()} sites")
    hg = pysam.FastaFile(a.genome_fasta)
    pos["base_tx"] = [hg.fetch(c, int(p), int(p) + 1).upper() for c, p in zip(pos["chrom"], pos["pos"])]
    pos.loc[pos["tx_strand"] == "-", "base_tx"] = pos.loc[pos["tx_strand"] == "-", "base_tx"].map(comp)

    # ------------------------------------------------------------------ (1) ancestral state of TE sites
    te = pos["set"] == "lnc_TE_site"
    if te.any():
        mp = pd.read_parquet(a.map_align, columns=["copy_id", "pos", "cons_pos"]).rename(columns={"cons_pos": "cons_pos_orig"})
        mp = mp[mp["copy_id"].isin(set(pos.loc[te, "copy_id"]))]
        pos = pos.merge(mp, on=["copy_id", "pos"], how="left")
        cons = load_consensus(a.consensus)
        te = pos["set"] == "lnc_TE_site"
        cb = []
        for sf, cpo, ts, xs in zip(pos["subfamily_orig"], pos["cons_pos_orig"], pos["te_strand"], pos["tx_strand"]):
            s = cons.get(str(sf), "")
            if not isinstance(cpo, (int, float)) or not np.isfinite(cpo) or not (0 < cpo <= len(s)):
                cb.append("N")
                continue
            b = s[int(cpo) - 1]
            cb.append(b if ts == xs else comp(b))
        pos["cons_tx"] = np.where(te, cb, "")
    else:
        pos["cons_tx"] = ""

    # ------------------------------------------------------------------ per site
    def agg_site(g):
        g = g.sort_values("dinuc_idx")
        af = (g["ac"] / g["an"]).where(g["callable"].astype(bool) & (g["an"] > 0))
        return pd.Series({"chrom": g["chrom"].iloc[0], "pos1": int(g["pos"].iloc[0]), "pos2": int(g["pos"].iloc[-1]),
                          "tx_strand": g["tx_strand"].iloc[0], "gene_id": g["gene_id"].iloc[0],
                          "site_type": g["site_type"].iloc[0], "dinuc": "".join(g["base_tx"]),
                          "cons_dinuc": "".join(g["cons_tx"].astype(str)),
                          "copy_id": g["copy_id"].iloc[0], "subfamily": g["subfamily"].iloc[0],
                          "subfamily_orig": g["subfamily_orig"].iloc[0], "rep_family": g["rep_family"].iloc[0],
                          "te_orientation": g["te_orientation"].iloc[0], "cons_pos": g["cons_pos"].iloc[0],
                          "max_af": float(np.nanmax(af)) if af.notna().any() else 0.0,
                          "polymorphic": bool((g["ac"] > 0).any())})
    S = pos.groupby(["set", "site_id"], sort=False).apply(agg_site, include_groups=False).reset_index()
    S["motif"] = S["site_type"].map(MOTIF)
    bad = S["dinuc"] != S["motif"]
    if bad.any():
        log(f"  WARNING: {int(bad.sum())} sites whose hg38 dinucleotide is not GT/AG (dropped)")
        S = S[~bad].copy()
    S["ancestral"] = np.where(S["set"] != "lnc_TE_site", "",
                              np.where(S["cons_dinuc"].str.contains("N"), "unknown",
                                       np.where(S["cons_dinuc"] == S["motif"], "preformed", "created")))
    S["segregating"] = S["polymorphic"]
    S["common_variant"] = S["max_af"] >= 0.01

    # ------------------------------------------------------------------ (2) other species
    from pyliftover import LiftOver
    copies = pd.read_csv(a.copies, sep="\t", usecols=["copy_id", "chrom", "start", "end"], dtype={"chrom": str})
    copies = copies.set_index("copy_id")
    names = []
    for item in a.species:
        name, chain, gpath = item.split("=")
        names.append(name)
        log(f"{name}: loading {chain}")
        lo = LiftOver(chain)
        gen = Genome(gpath)
        st, cpres = [], []
        for r in S.itertuples(index=False):
            m1, m2 = lift_one(lo, r.chrom, r.pos1), lift_one(lo, r.chrom, r.pos2)
            if m1 is None or m2 is None or m1[0] != m2[0] or m1[2] != m2[2] or abs(m1[1] - m2[1]) != 1:
                st.append("unaligned")
            else:
                flip = (m1[2] == "-") != (r.tx_strand == "-")
                b1, b2 = gen.base(m1[0], m1[1]), gen.base(m2[0], m2[1])
                if flip:
                    b1, b2 = comp(b1), comp(b2)
                st.append("present" if b1 + b2 == r.motif else ("absent" if "N" not in b1 + b2 else "unaligned"))
            if r.set == "lnc_TE_site" and isinstance(r.copy_id, str) and r.copy_id in copies.index:
                c = copies.loc[r.copy_id]
                L = int(c["end"]) - int(c["start"]) - 2 * a.copy_trim
                if L < 20:
                    cpres.append("short")
                    continue
                e1 = lift_one(lo, c["chrom"], int(c["start"]) + a.copy_trim)
                e2 = lift_one(lo, c["chrom"], int(c["end"]) - 1 - a.copy_trim)
                ok = (e1 is not None and e2 is not None and e1[0] == e2[0] and e1[2] == e2[2]
                      and 0.5 * L <= abs(e2[1] - e1[1]) <= 2 * L)
                cpres.append("present" if ok else "absent")
            else:
                cpres.append("")
        S[f"site_{name}"] = st
        S[f"copy_{name}"] = cpres
        log(f"  {name}: " + ", ".join(f"{k} {v}" for k, v in S[f"site_{name}"].value_counts().items()))
        del lo

    # ------------------------------------------------------------------ (3) age class
    if names:
        pres = np.stack([S[f"site_{n}"].to_numpy() == "present" for n in names], axis=1)
        alig = np.stack([S[f"site_{n}"].to_numpy() != "unaligned" for n in names], axis=1)
        deepest = np.where(pres.any(axis=1), pres.shape[1] - 1 - np.argmax(pres[:, ::-1], axis=1), -1)
        age = np.where(deepest >= 0, np.array([f"shared with {n}" for n in names], dtype=object)[np.maximum(deepest, 0)],
                       "unresolved")
        if len(names) >= 2:
            hs = (deepest < 0) & alig[:, 0] & alig[:, 1]
            age = np.where(hs, "human-specific", age)
        S["age_class"] = age
        S["age_rank"] = np.where(deepest >= 0, deepest + 1, np.where(age == "human-specific", 0, -1))
        # TE-derived sites in copies shared with the closest species but absent there: created after the split
        n0 = names[0]
        S[f"created_after_{n0}_split"] = np.where(S["set"] == "lnc_TE_site",
                                                 (S[f"copy_{n0}"] == "present") & (S[f"site_{n0}"] == "absent"), False)
    S.to_csv(f"{a.out}_sites.tsv.gz", sep="\t", index=False)

    # ------------------------------------------------------------------ summaries
    rows = []
    def groups():
        for st, s in S.groupby("set"):
            yield st, "all", s
            for col in ("site_type",):
                for v, x in s.groupby(col):
                    yield st, f"{col}={v}", x
            if st == "lnc_TE_site":
                for col in ("rep_family", "te_orientation", "ancestral"):
                    for v, x in s.groupby(col):
                        yield st, f"{col}={v}", x
    for st, sub, s in groups():
        base = {"set": st, "subset": sub, "sites": len(s)}
        cats = [("segregating", s["segregating"]), ("common_variant", s["common_variant"])]
        if st == "lnc_TE_site":
            cats.append(("preformed", s["ancestral"] == "preformed"))
        for n in names:
            al = s[f"site_{n}"] != "unaligned"
            k, m = int((s[f"site_{n}"] == "present").sum()), int(al.sum())
            p, lo_, hi_ = wilson(k, m)
            rows.append({**base, "measure": f"present_in_{n} (of aligned)", "n": m, "share": p, "lo": lo_, "hi": hi_})
            rows.append({**base, "measure": f"aligned_in_{n}", "n": len(s), "share": al.mean(), "lo": np.nan, "hi": np.nan})
        if names:
            for cl, x in s["age_class"].value_counts().items():
                p, lo_, hi_ = wilson(int(x), len(s))
                rows.append({**base, "measure": f"age={cl}", "n": len(s), "share": p, "lo": lo_, "hi": hi_})
            if st == "lnc_TE_site":
                c0 = f"created_after_{names[0]}_split"
                sh = s[s[f"copy_{names[0]}"] == "present"]
                p, lo_, hi_ = wilson(int(sh[c0].sum()), len(sh))
                rows.append({**base, "measure": f"{c0} (of sites in shared copies)", "n": len(sh), "share": p,
                             "lo": lo_, "hi": hi_})
        for name, m in cats:
            p, lo_, hi_ = wilson(int(m.sum()), len(s))
            rows.append({**base, "measure": name, "n": len(s), "share": p, "lo": lo_, "hi": hi_})
    T = pd.DataFrame(rows)
    T.to_csv(f"{a.out}_summary.tsv", sep="\t", index=False)
    pd.set_option("display.width", 220)
    log(T[T["subset"] == "all"].round(4).to_string(index=False))


if __name__ == "__main__":
    main()
