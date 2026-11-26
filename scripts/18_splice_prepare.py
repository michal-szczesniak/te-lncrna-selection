#!/usr/bin/env python3
"""Prepare the splice-site constraint test and the LTR-promoter table.

Splice sites: canonical GT-AG introns of all annotated transcripts of the kept lncRNA genes and of kept
protein-coding genes (step-1 gene tables). For each intron the two donor bases (GT) and the two acceptor
bases (AG) are taken. Sets written to positions.tsv.gz:
  lnc_TE_site       dinucleotides lying inside a TE copy mapped in step 3 (exact consensus position)
  lnc_TE_ctrl       matched controls: the same subfamily, consensus position and base (in TE orientation)
                    in mapped intronic or intergenic copies that host no splice site (--ctrl-per-site each)
  lnc_nonTE_site    lncRNA dinucleotides outside any repeat (random --max-introns introns)
  lnc_nonTE_flank   control bases --flank nt deeper in the same introns
  pc_site, pc_flank the same for protein-coding genes (random --max-introns introns)
Positions of TE sites and controls are inside the step-2 gnomAD regions (already fetched); the others are
written to fetch_regions.bed for a new gnomAD extraction.

LTR promoters: transcription starts of lncRNA and protein-coding transcripts inside LTR elements, by
RepeatMasker subfamily, with obs/exp against random exonic positions of the same genes (ltr_tss.tsv).
"""
import argparse
import os

import numpy as np
import pandas as pd
import pysam

from ridl_utils import index_by_chrom, log, merge_intervals, open_text, read_gtf

COMP = str.maketrans("ACGTN", "TGCAN")


def read_rmsk(path, chroms):
    rows = []
    with open_text(path) as fh:
        for line in fh:
            f = line.split()
            if len(f) < 15 or not f[0].isdigit() or f[4] not in chroms:
                continue
            cf = f[10].split("/")
            rows.append((f[4], int(f[5]) - 1, int(f[6]), "+" if f[8] == "+" else "-", f[9], cf[0], cf[-1]))
    r = pd.DataFrame(rows, columns=["chrom", "start", "end", "strand", "name", "rclass", "rfamily"])
    return {c: g.sort_values("start").reset_index(drop=True) for c, g in r.groupby("chrom")}


def hit_rmsk(chrom, P, rm):
    """index of the repeat record covering each position (or -1)."""
    hit = np.full(len(P), -1)
    if chrom not in rm:
        return hit
    S, E = rm[chrom]["start"].to_numpy(), rm[chrom]["end"].to_numpy()
    idx = np.searchsorted(S, P, side="right") - 1
    for k in range(6):
        cand = idx - k
        ok = (hit < 0) & (cand >= 0)
        ok[ok] &= E[cand[ok]] > P[ok]
        hit[ok] = cand[ok]
    return hit


def introns_of(tx, ex):
    e = ex.sort_values(["transcript_id", "start"])
    nxt = e.groupby("transcript_id")["start"].shift(-1)
    m = nxt.notna()
    it = pd.DataFrame({"gene_id": e.loc[m, "gene_id"], "chrom": e.loc[m, "chrom"], "strand": e.loc[m, "strand"],
                       "start": e.loc[m, "end"].to_numpy(), "end": nxt[m].astype(np.int64).to_numpy()})
    it = it[it["end"] - it["start"] >= 30].drop_duplicates(["chrom", "start", "end", "strand"])
    return it.reset_index(drop=True)


def site_positions(it, fa, flank):
    """Rows: one per base of donor/acceptor dinucleotides (+ flank bases); canonical GT-AG only."""
    rows = []
    for k, r in enumerate(it.itertuples(index=False)):
        s, e = int(r.start), int(r.end)
        if r.strand == "+":
            don, acc = [s, s + 1], [e - 2, e - 1]
            fdon, facc = [s + flank, s + flank + 1], [e - flank - 2, e - flank - 1]
        else:
            don, acc = [e - 1, e - 2], [s + 1, s]
            fdon, facc = [e - flank - 1, e - flank - 2], [s + flank + 1, s + flank]
        seq = lambda ps: "".join(fa.fetch(r.chrom, p, p + 1).upper() for p in ps)
        dseq, aseq = seq(don), seq(acc)
        if r.strand == "-":
            dseq, aseq = dseq.translate(COMP), aseq.translate(COMP)
        if dseq != "GT" or aseq != "AG":
            continue
        sid = f"{r.chrom}:{s}-{e}{r.strand}"
        long_enough = (e - s) >= 2 * flank + 10
        for typ, ps, fps in (("donor", don, fdon), ("acceptor", acc, facc)):
            for j, p in enumerate(ps):
                rows.append((r.chrom, p, sid + typ[0], typ, j + 1, "site", r.strand, r.gene_id))
            if long_enough:
                for j, p in enumerate(fps):
                    rows.append((r.chrom, p, sid + typ[0], typ, j + 1, "flank", r.strand, r.gene_id))
    return pd.DataFrame(rows, columns=["chrom", "pos", "site_id", "site_type", "dinuc_idx", "role", "tx_strand",
                                       "gene_id"])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--lnc-genes", required=True)
    ap.add_argument("--pc-genes", required=True)
    ap.add_argument("--rmsk-out", required=True)
    ap.add_argument("--genome-fasta", required=True)
    ap.add_argument("--copies", required=True, help="step-3b copies_grouped.tsv.gz")
    ap.add_argument("--map", required=True, help="step-3b map_grouped.parquet")
    ap.add_argument("--te-regions", required=True, help="step-2 gnomad_regions.bed (already fetched)")
    ap.add_argument("--max-introns", type=int, default=10000)
    ap.add_argument("--ctrl-per-site", type=int, default=5)
    ap.add_argument("--flank", type=int, default=60)
    ap.add_argument("--bg-per-gene", type=int, default=50)
    ap.add_argument("--chroms", default=",".join(f"chr{i}" for i in range(1, 23)))
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out-dir", required=True)
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    rng = np.random.default_rng(a.seed)
    chroms = a.chroms.split(",")
    fa = pysam.FastaFile(a.genome_fasta)
    genes, tx, ex = read_gtf(a.gtf, chroms)
    rm = read_rmsk(a.rmsk_out, set(chroms))
    kept = {}
    for name, path in (("lnc", a.lnc_genes), ("pc", a.pc_genes)):
        g = pd.read_csv(path, sep="\t")
        kept[name] = set(g.loc[g["drop_reason"].fillna("") == "", "gene_id"])

    # ------------------------------------------------------------- splice sites
    pos_tables = []
    copies = pd.read_csv(a.copies, sep="\t", usecols=["copy_id", "chrom", "start", "end", "te_strand", "subfamily",
                                                      "subfamily_orig", "rep_family", "compartment"])
    mp = pd.read_parquet(a.map).merge(copies[["copy_id", "chrom"]], on="copy_id")
    mp_idx = mp.set_index(["chrom", "pos"])
    for name in ("lnc", "pc"):
        t = tx[tx["gene_id"].isin(kept[name])]
        e = ex[ex["gene_id"].isin(kept[name])]
        it = introns_of(t, e)
        log(f"{name}: introns {len(it)}")
        if name == "pc" and len(it) > 2 * a.max_introns:
            it = it.iloc[np.sort(rng.choice(len(it), 2 * a.max_introns, replace=False))].reset_index(drop=True)
        sp = site_positions(it, fa, a.flank)
        log(f"{name}: canonical GT-AG introns {sp['site_id'].str[:-1].nunique()}")
        # repeat status of the dinucleotide bases
        sp["in_rep"] = False
        for c, g in sp.groupby("chrom"):
            sp.loc[g.index, "in_rep"] = hit_rmsk(c, g["pos"].to_numpy(), rm) >= 0
        if name == "pc":
            ids = sp["site_id"].str[:-1].unique()
            pick = set(rng.choice(ids, min(a.max_introns, len(ids)), replace=False))
            s = sp[sp["site_id"].str[:-1].isin(pick)].copy()
            s["set"] = np.where(s["role"] == "site", "pc_site", "pc_flank")
            pos_tables.append(s)
            continue
        # lncRNA: TE sites (both dinucleotide bases in a mapped copy) and non-TE sites (no repeat at all)
        site_rows = sp[sp["role"] == "site"]
        key = pd.MultiIndex.from_frame(site_rows[["chrom", "pos"]])
        inmap = key.isin(mp_idx.index)
        site_rows = site_rows.assign(inmap=inmap)
        per_site = site_rows.groupby("site_id").agg(n_map=("inmap", "sum"), n_rep=("in_rep", "sum"))
        te_ids = set(per_site.index[per_site["n_map"] == 2])
        non_ids = per_site.index[per_site["n_rep"] == 0]
        non_ids = set(rng.choice(non_ids, min(a.max_introns * 2, len(non_ids)), replace=False))
        te = site_rows[site_rows["site_id"].isin(te_ids)].copy()
        te = te.merge(mp[["chrom", "pos", "copy_id", "cons_pos"]], on=["chrom", "pos"], how="left")
        te = te.drop_duplicates(["chrom", "pos", "site_id"]).merge(
            copies[["copy_id", "te_strand", "subfamily", "subfamily_orig", "rep_family"]], on="copy_id")
        te["set"] = "lnc_TE_site"
        te["te_orientation"] = np.where(te["te_strand"] == te["tx_strand"], "sense", "antisense")
        non = sp[sp["site_id"].isin(non_ids)].copy()
        non["set"] = np.where(non["role"] == "site", "lnc_nonTE_site", "lnc_nonTE_flank")
        log(f"lnc: TE-derived sites (dinucleotide inside a mapped copy) {len(te_ids)}, "
            f"non-TE sites kept {len(non_ids)}, sites in other repeats ignored "
            f"{int(((per_site['n_rep'] > 0) & (per_site['n_map'] < 2)).sum())}")
        pos_tables += [te, non]

    pos = pd.concat(pos_tables, ignore_index=True)
    # base in TE orientation for TE sites -> stratum
    te = pos[pos["set"] == "lnc_TE_site"].copy()
    te["base_g"] = [fa.fetch(c, int(p), int(p) + 1).upper() for c, p in zip(te["chrom"], te["pos"])]
    te["base_te"] = np.where(te["te_strand"] == "+", te["base_g"], te["base_g"].str.translate(COMP))
    te["stratum"] = te["subfamily_orig"] + "|" + te["cons_pos"].astype(int).astype(str) + "|" + te["base_te"]

    # ------------------------------------------------------------- matched controls
    site_keys = pd.MultiIndex.from_frame(pos[["chrom", "pos"]])
    hosting = mp.set_index(["chrom", "pos"]).loc[lambda d: d.index.isin(site_keys), "copy_id"].unique()
    pool_copies = copies[copies["compartment"].isin(["intron", "intergenic"]) & ~copies["copy_id"].isin(hosting)]
    strata = te[["subfamily_orig", "cons_pos"]].drop_duplicates()
    pool = mp[mp["copy_id"].isin(set(pool_copies["copy_id"]))].merge(
        pool_copies[["copy_id", "te_strand", "subfamily", "subfamily_orig", "rep_family"]], on="copy_id")
    pool = pool.merge(strata, on=["subfamily_orig", "cons_pos"])
    pool["base_g"] = [fa.fetch(c, int(p), int(p) + 1).upper() for c, p in zip(pool["chrom"], pool["pos"])]
    pool["base_te"] = np.where(pool["te_strand"] == "+", pool["base_g"], pool["base_g"].str.translate(COMP))
    pool["stratum"] = pool["subfamily_orig"] + "|" + pool["cons_pos"].astype(int).astype(str) + "|" + pool["base_te"]
    need = te.groupby("stratum").size() * a.ctrl_per_site
    ctrl = []
    for st, g in pool[pool["stratum"].isin(need.index)].groupby("stratum"):
        k = min(int(need[st]), len(g))
        ctrl.append(g.iloc[rng.choice(len(g), k, replace=False)])
    ctrl = pd.concat(ctrl, ignore_index=True) if ctrl else pool.iloc[:0]
    ctrl["set"] = "lnc_TE_ctrl"
    ctrl["site_id"] = ctrl["stratum"]
    log(f"TE sites: {te['site_id'].nunique()} sites, {len(te)} bases in {te['stratum'].nunique()} strata; "
        f"controls: {len(ctrl)} bases in {ctrl['stratum'].nunique()} strata")
    out = pd.concat([pos[pos["set"] != "lnc_TE_site"], te, ctrl], ignore_index=True)
    cols = ["chrom", "pos", "set", "site_id", "site_type", "dinuc_idx", "tx_strand", "gene_id", "copy_id", "te_strand",
            "te_orientation", "subfamily", "subfamily_orig", "rep_family", "cons_pos", "stratum"]
    for c in cols:
        if c not in out.columns:
            out[c] = np.nan
    out[cols].to_csv(os.path.join(a.out_dir, "positions.tsv.gz"), sep="\t", index=False)

    # ------------------------------------------------------------- regions to fetch from gnomAD
    tr = pd.read_csv(a.te_regions, sep="\t", header=None, usecols=[0, 1, 2], names=["chrom", "start", "end"])
    tidx = index_by_chrom(merge_intervals(tr))
    out["in_te_regions"] = False
    for c, g in out.groupby("chrom"):
        if c in tidx:
            S, E = tidx[c]
            i = np.searchsorted(E, g["pos"].to_numpy(), side="right")
            ok = i < len(S)
            inside = np.zeros(len(g), dtype=bool)
            inside[ok] = S[i[ok]] <= g["pos"].to_numpy()[ok]
            out.loc[g.index, "in_te_regions"] = inside
    f = out[~out["in_te_regions"]][["chrom", "pos"]].drop_duplicates()
    bed = merge_intervals(pd.DataFrame({"chrom": f["chrom"], "start": f["pos"], "end": f["pos"] + 1}))
    bed.sort_values(["chrom", "start"]).to_csv(os.path.join(a.out_dir, "fetch_regions.bed"), sep="\t",
                                                index=False, header=False)
    summ = out.groupby(["set", "in_te_regions"]).size().rename("bases").reset_index()
    summ.to_csv(os.path.join(a.out_dir, "prepare_summary.tsv"), sep="\t", index=False)
    log(summ.to_string(index=False))
    log(f"regions to fetch from gnomAD: {len(bed)}")

    # ------------------------------------------------------------- LTR promoters
    rows = []
    for name in ("lnc", "pc"):
        t = tx[tx["gene_id"].isin(kept[name])]
        tss = pd.DataFrame({"gene_id": t["gene_id"], "chrom": t["chrom"],
                            "pos": np.where(t["strand"] == "+", t["start"], t["end"] - 1)}).drop_duplicates()
        em = ex[ex["gene_id"].isin(kept[name])].drop_duplicates(["gene_id", "start", "end"])
        bg = []
        for gid, g in em.groupby("gene_id"):
            L = (g["end"] - g["start"]).to_numpy()
            k = rng.choice(len(g), a.bg_per_gene, p=L / L.sum())
            bg.append(pd.DataFrame({"chrom": g["chrom"].iloc[0],
                                    "pos": g["start"].to_numpy()[k] + (rng.random(a.bg_per_gene) * L[k]).astype(np.int64)}))
        bg = pd.concat(bg, ignore_index=True)
        for what, d in (("TSS", tss), ("bg_exon", bg)):
            names = []
            for c, g in d.groupby("chrom"):
                h = hit_rmsk(c, g["pos"].to_numpy(), rm)
                r = rm.get(c)
                for x in h:
                    names.append((r.at[x, "rclass"], r.at[x, "rfamily"], r.at[x, "name"]) if x >= 0 else ("", "", ""))
            nm = pd.DataFrame(names, columns=["rclass", "rfamily", "name"])
            nm = nm[nm["rclass"] == "LTR"]
            cnt = nm.groupby(["rfamily", "name"]).size().rename("n").reset_index()
            cnt["share"] = cnt["n"] / len(d)
            cnt["gene_set"], cnt["what"] = name, what
            rows.append(cnt)
    lt = pd.concat(rows, ignore_index=True)
    w = lt.pivot_table(index=["rfamily", "name"], columns=["gene_set", "what"], values="share", fill_value=0)
    w.columns = [f"{g}_{x}" for g, x in w.columns]
    n = lt[(lt.gene_set == "lnc") & (lt.what == "TSS")].set_index(["rfamily", "name"])["n"]
    w["lnc_TSS_n"] = n.reindex(w.index).fillna(0).astype(int)
    if "lnc_bg_exon" in w:
        w["lnc_TSS_obs_exp"] = w["lnc_TSS"] / w["lnc_bg_exon"].replace(0, np.nan)
    if "pc_TSS" in w:
        w["lnc_vs_pc_TSS"] = w["lnc_TSS"] / w["pc_TSS"].replace(0, np.nan)
    w = w.sort_values("lnc_TSS_n", ascending=False).reset_index()
    w.to_csv(os.path.join(a.out_dir, "ltr_tss.tsv"), sep="\t", index=False)
    fam = lt[lt.what == "TSS"].groupby(["gene_set", "rfamily"])["share"].sum().unstack(0)
    log("LTR families at TSS (share of all TSS):\n" + fam.round(4).to_string())
    pd.set_option("display.width", 200)
    log("top LTR subfamilies at lncRNA TSS:\n" + w.head(15).round(4).to_string(index=False))


if __name__ == "__main__":
    main()
