#!/usr/bin/env python3
"""Revision step 18b. Alternative controls for TE-derived splice sites (location-matched).

Step 18 compares each TE-derived splice dinucleotide with the same consensus position and base (stratum
subfamily_orig|cons_pos|base) in unused intronic or intergenic copies anywhere in the genome. Because the
paper shows that genomic location alone changes variant density, this script draws controls closer in
location:
  --pool intron      unused copies whose (dominant) compartment is intron
  --pool samegene    unused copies (any compartment except intergenic) in the SAME gene as the site; each
                     site base is matched only with bases of its own gene (stratum x gene)
Controls are positions in step-2/3 copies, so their gnomAD data are already in the step-2 subsets.

Input: positions.tsv.gz of step 18 (TE sites and their strata), step-3b map and copies, genome FASTA.
Output: <out-dir>/positions.tsv.gz with the TE sites and the new controls (set lnc_TE_ctrl), ready for
19_splice_analyze.py, and a summary of how many sites found controls.
"""
import argparse
import os

import numpy as np
import pandas as pd
import pysam

from ridl_utils import log

COMP = str.maketrans("ACGTN", "TGCAN")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--positions", required=True, help="step-18 positions.tsv.gz")
    ap.add_argument("--copies", required=True)
    ap.add_argument("--map", required=True)
    ap.add_argument("--genome-fasta", required=True)
    ap.add_argument("--pool", choices=["intron", "samegene"], required=True)
    ap.add_argument("--ctrl-per-site", type=int, default=20)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out-dir", required=True)
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    rng = np.random.default_rng(a.seed)
    fa = pysam.FastaFile(a.genome_fasta)
    pos = pd.read_csv(a.positions, sep="\t", low_memory=False)
    te = pos[pos["set"] == "lnc_TE_site"].copy()
    copies = pd.read_csv(a.copies, sep="\t", usecols=["copy_id", "chrom", "te_strand", "subfamily", "subfamily_orig",
                                                      "rep_family", "compartment", "gene_id"], dtype={"chrom": str})
    mp = pd.read_parquet(a.map, columns=["copy_id", "pos", "cons_pos"])
    hosting = set(pos.loc[pos["set"] == "lnc_TE_site", "copy_id"].dropna())
    if a.pool == "intron":
        pc = copies[(copies["compartment"] == "intron") & ~copies["copy_id"].isin(hosting)]
    else:
        pc = copies[(copies["compartment"] != "intergenic") & ~copies["copy_id"].isin(hosting)
                    & copies["gene_id"].isin(set(te["gene_id"]))]
    strata = te[["subfamily_orig", "cons_pos"]].drop_duplicates()
    pool = mp[mp["copy_id"].isin(set(pc["copy_id"]))].merge(pc, on="copy_id")
    pool = pool.merge(strata, on=["subfamily_orig", "cons_pos"])
    pool["base_g"] = [fa.fetch(c, int(p), int(p) + 1).upper() for c, p in zip(pool["chrom"], pool["pos"])]
    pool["base_te"] = np.where(pool["te_strand"] == "+", pool["base_g"], pool["base_g"].str.translate(COMP))
    pool["stratum"] = pool["subfamily_orig"] + "|" + pool["cons_pos"].astype(int).astype(str) + "|" + pool["base_te"]
    ctrl = []
    if a.pool == "intron":
        need = te.groupby("stratum").size() * a.ctrl_per_site
        for st, g in pool[pool["stratum"].isin(need.index)].groupby("stratum"):
            k = min(int(need[st]), len(g))
            ctrl.append(g.iloc[rng.choice(len(g), k, replace=False)])
        keep_te = te
    else:
        te["_key"] = te["stratum"] + "@" + te["gene_id"].astype(str)
        pool["_key"] = pool["stratum"] + "@" + pool["gene_id"].astype(str)
        need = te.groupby("_key").size() * a.ctrl_per_site
        for key, g in pool[pool["_key"].isin(need.index)].groupby("_key"):
            k = min(int(need[key]), len(g))
            ctrl.append(g.iloc[rng.choice(len(g), k, replace=False)])
        found = set(pool["_key"])
        keep_te = te[te["_key"].isin(found)].copy()
        # a site enters only if both of its bases found a same-gene control
        ok = keep_te.groupby("site_id")["_key"].transform("size") == 2
        keep_te = keep_te[ok].drop(columns="_key")
        ctrl = [c.drop(columns="_key") for c in ctrl]
    ctrl = pd.concat(ctrl, ignore_index=True) if ctrl else pool.iloc[:0]
    ctrl["set"] = "lnc_TE_ctrl"
    ctrl["site_id"] = ctrl["stratum"]
    for c in ("site_type", "dinuc_idx", "tx_strand", "te_orientation"):
        ctrl[c] = np.nan
    cols = ["chrom", "pos", "set", "site_id", "site_type", "dinuc_idx", "tx_strand", "gene_id", "copy_id", "te_strand",
            "te_orientation", "subfamily", "subfamily_orig", "rep_family", "cons_pos", "stratum"]
    out = pd.concat([keep_te[cols], ctrl[cols]], ignore_index=True)
    out.to_csv(os.path.join(a.out_dir, "positions.tsv.gz"), sep="\t", index=False)
    msg = (f"pool {a.pool}: TE sites kept {keep_te['site_id'].nunique()} of {te['site_id'].nunique()}, "
           f"control bases {len(ctrl)} in {ctrl['stratum'].nunique()} strata")
    log(msg)
    with open(os.path.join(a.out_dir, "summary.txt"), "w") as fh:
        fh.write(msg + "\n")


if __name__ == "__main__":
    main()
