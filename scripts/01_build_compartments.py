#!/usr/bin/env python3
"""Step 1. Build genomic compartments for lncRNA genes.

Outputs (in --out-dir):
  compartments.tsv.gz   chrom,start,end,compartment(exon|intron|intergenic),gene_id,gene_name,
                        strand,exon_class(mono|multi),max_tpm   -- non-overlapping intervals
  splice_sites.tsv.gz   chrom,pos,gene_id   (exon/intron boundaries of kept lncRNA transcripts)
  ref_tx.tsv.gz         reference transcript exons per kept gene, with 5'->3' transcript offsets
  lnc_genes.tsv.gz      kept genes and why others were dropped (summary)
"""
import argparse
import os

import numpy as np
import pandas as pd

from ridl_utils import (index_by_chrom, log, merge_intervals, overlaps_any, parse_chroms,
                        points_to_intervals, read_bed3, read_gtf, subtract)


def splice_and_ends(tx, ex):
    """Return splice-site points (exon boundaries internal to transcripts) + TSS/TES points."""
    ex = ex.sort_values(["transcript_id", "start"])
    n_ex = ex.groupby("transcript_id")["start"].transform("size")
    first = ~ex["transcript_id"].duplicated(keep="first")
    last = ~ex["transcript_id"].duplicated(keep="last")
    multi = n_ex > 1
    # internal boundaries: start of non-first exons, end of non-last exons
    s_pts = ex.loc[multi & ~first, ["chrom", "start", "gene_id"]].rename(columns={"start": "pos"})
    e_pts = ex.loc[multi & ~last, ["chrom", "end", "gene_id"]].rename(columns={"end": "pos"})
    splice = pd.concat([s_pts, e_pts], ignore_index=True).drop_duplicates()
    tss = np.where(tx["strand"] == "+", tx["start"], tx["end"])
    tes = np.where(tx["strand"] == "+", tx["end"], tx["start"])
    ends = pd.DataFrame({"chrom": tx["chrom"].to_numpy(), "tss": tss, "tes": tes,
                         "gene_id": tx["gene_id"].to_numpy()})
    return splice, ends


def read_gtex_max_tpm(path):
    df = pd.read_csv(path, sep="\t", skiprows=2)
    name_col = "Name" if "Name" in df.columns else df.columns[0]
    num = df.select_dtypes(include=[np.number])
    out = pd.DataFrame({"gid": df[name_col].astype(str).str.replace(r"\.\d+(_PAR_Y)?$", "", regex=True),
                        "max_tpm": num.max(axis=1)})
    return out.groupby("gid")["max_tpm"].max()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--chrom-sizes", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--chroms", default=None, help="comma list; default chr1-22,chrX")
    ap.add_argument("--lnc-types", default="lncRNA",
                    help="gene types to build compartments for (e.g. protein_coding for a positive control)")
    ap.add_argument("--skip-overlap-filter", action="store_true",
                    help="do not drop genes whose exons overlap protein-coding/pseudogene exons (needed when "
                         "--lnc-types protein_coding)")
    ap.add_argument("--gtex-gct", default=None, help="GTEx gene median TPM .gct(.gz)")
    ap.add_argument("--min-tpm", type=float, default=1.0)
    ap.add_argument("--keep-unmatched-gtex", action="store_true",
                    help="keep genes absent from the GTEx table (default: drop)")
    ap.add_argument("--mask-bed", action="append", default=[],
                    help="BED of regions removed from all compartments (e.g. segmental dups)")
    ap.add_argument("--exon-splice-pad", type=int, default=10)
    ap.add_argument("--intron-splice-pad", type=int, default=50)
    ap.add_argument("--exon-tss-pad", type=int, default=500)
    ap.add_argument("--intron-tss-pad", type=int, default=1000)
    ap.add_argument("--tes-pad", type=int, default=100)
    ap.add_argument("--intergenic-gene-pad", type=int, default=10000)
    ap.add_argument("--other-genes", choices=["spans", "exons"], default="spans",
                    help="remove from a gene's compartments the whole spans of other genes (default) or only "
                         "their exons; with 'exons', intronic sequence shared by two kept genes is also removed")
    ap.add_argument("--keep-other-gene-spans", action="store_true",
                    help="do NOT remove regions covered by other genes' spans (default: remove)")
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    chroms = parse_chroms(a.chroms)

    log("reading GTF ...")
    genes, tx, ex = read_gtf(a.gtf, chroms)
    lnc_types = set(a.lnc_types.split(","))
    lnc = genes[genes["gene_type"].isin(lnc_types)].copy()
    log(f"  genes: {len(genes)}, lncRNA genes: {len(lnc)}")

    # --- drop lncRNAs whose exons overlap protein-coding or pseudogene exons (any strand)
    bad_types = genes["gene_type"].eq("protein_coding") | genes["gene_type"].str.contains("pseudogene", na=False)
    bad_ex = ex[ex["gene_id"].isin(genes.loc[bad_types, "gene_id"])]
    bad_idx = index_by_chrom(merge_intervals(bad_ex[["chrom", "start", "end"]]))
    lnc_ex = ex[ex["gene_id"].isin(lnc["gene_id"])].reset_index(drop=True)
    hit = overlaps_any(lnc_ex[["chrom", "start", "end"]], bad_idx)
    drop_pc = set() if a.skip_overlap_filter else set(lnc_ex.loc[hit, "gene_id"])
    lnc["drop_reason"] = np.where(lnc["gene_id"].isin(drop_pc), "exon_overlaps_pc_or_pseudo", "")

    # --- expression filter
    lnc["max_tpm"] = np.nan
    if a.gtex_gct:
        tpm = read_gtex_max_tpm(a.gtex_gct)
        gid = lnc["gene_id"].str.replace(r"\.\d+(_PAR_Y)?$", "", regex=True)
        lnc["max_tpm"] = gid.map(tpm).to_numpy()
        unmatched = lnc["max_tpm"].isna()
        low = lnc["max_tpm"] < a.min_tpm
        lnc.loc[(lnc["drop_reason"] == "") & low, "drop_reason"] = "low_expression"
        if not a.keep_unmatched_gtex:
            lnc.loc[(lnc["drop_reason"] == "") & unmatched, "drop_reason"] = "not_in_gtex"
    kept = lnc[lnc["drop_reason"] == ""].copy()
    log(f"  kept lncRNA genes: {len(kept)}  (dropped: {lnc['drop_reason'].value_counts().to_dict()})")
    lnc.to_csv(os.path.join(a.out_dir, "lnc_genes.tsv.gz"), sep="\t", index=False)

    # --- exon class (mono/multi) per kept gene
    kept_ex = ex[ex["gene_id"].isin(kept["gene_id"])]
    n_ex_tx = kept_ex.groupby(["gene_id", "transcript_id"]).size().groupby("gene_id").max()
    kept["exon_class"] = np.where(kept["gene_id"].map(n_ex_tx) > 1, "multi", "mono")

    # --- points
    splice_all, ends_all = splice_and_ends(tx, ex)
    splice_kept = splice_all[splice_all["gene_id"].isin(kept["gene_id"])]

    masks = read_bed3(a.mask_bed)
    masks = masks[masks["chrom"].isin(chroms)]

    def merged_index(*dfs):
        d = pd.concat([x[["chrom", "start", "end"]] for x in dfs if len(x)], ignore_index=True)
        return index_by_chrom(merge_intervals(d)) if len(d) else {}

    # other-gene spans: for each kept gene, spans of all *other* genes. Implemented by
    # counting coverage: positions covered by >=2 gene spans, where one is the kept gene,
    # are removed. Approximation: remove spans of all genes except the gene itself.
    all_spans = genes[["chrom", "start", "end", "gene_id"]]

    # --- exonic compartment
    log("building exonic compartment ...")
    ex_union = merge_intervals(kept_ex[["chrom", "start", "end", "gene_id"]], by=("gene_id", "chrom"))
    ex_union = ex_union.merge(kept[["gene_id", "gene_name", "strand", "exon_class", "max_tpm"]], on="gene_id")
    ex_union = ex_union[["chrom", "start", "end", "gene_id", "gene_name", "strand", "exon_class", "max_tpm"]]
    rm_exon = merged_index(
        points_to_intervals(splice_all["chrom"], splice_all["pos"], a.exon_splice_pad),
        points_to_intervals(ends_all["chrom"], ends_all["tss"], a.exon_tss_pad),
        points_to_intervals(ends_all["chrom"], ends_all["tes"], a.tes_pad),
        masks)
    exon_c = subtract(ex_union, rm_exon)

    # --- intronic compartment
    log("building intronic compartment ...")
    spans = kept[["chrom", "start", "end", "gene_id", "gene_name", "strand", "exon_class", "max_tpm"]]
    rm_intron = merged_index(
        ex[["chrom", "start", "end"]],
        points_to_intervals(splice_all["chrom"], splice_all["pos"], a.intron_splice_pad),
        points_to_intervals(ends_all["chrom"], ends_all["tss"], a.intron_tss_pad),
        masks)
    intron_c = subtract(spans, rm_intron)

    # --- remove regions covered by other genes' spans (or only their exons)
    def remove_other(comp, others_all):
        pieces = []
        for chrom, g in comp.groupby("chrom"):
            others = others_all[others_all["chrom"] == chrom]
            for gid, gg in g.groupby("gene_id"):
                o = others[others["gene_id"] != gid]
                lo, hi = gg["start"].min(), gg["end"].max()
                o = o[(o["end"] > lo) & (o["start"] < hi)]
                if len(o):
                    gg = subtract(gg, index_by_chrom(merge_intervals(o[["chrom", "start", "end"]])))
                pieces.append(gg)
        return pd.concat(pieces, ignore_index=True) if pieces else comp.iloc[0:0]

    if not a.keep_other_gene_spans and a.other_genes == "exons":
        log("removing other genes' exons; intronic sequence shared by kept genes ...")
        exon_c = remove_other(exon_c, ex[["chrom", "start", "end", "gene_id"]])
        intron_c = remove_other(intron_c, intron_c[["chrom", "start", "end", "gene_id"]])
    elif not a.keep_other_gene_spans:
        log("removing regions covered by other genes ...")
        parts = []
        for comp in (exon_c, intron_c):
            pieces = []
            for chrom, g in comp.groupby("chrom"):
                others = all_spans[all_spans["chrom"] == chrom]
                for gid, gg in g.groupby("gene_id"):
                    o = others[others["gene_id"] != gid]
                    lo, hi = gg["start"].min(), gg["end"].max()
                    o = o[(o["end"] > lo) & (o["start"] < hi)]
                    if len(o):
                        gg = subtract(gg, index_by_chrom(merge_intervals(o[["chrom", "start", "end"]])))
                    pieces.append(gg)
            parts.append(pd.concat(pieces, ignore_index=True) if pieces else comp.iloc[0:0])
        exon_c, intron_c = parts

    # --- intergenic compartment
    log("building intergenic compartment ...")
    sizes = pd.read_csv(a.chrom_sizes, sep="\t", header=None, names=["chrom", "size"], usecols=[0, 1])
    sizes = sizes[sizes["chrom"].isin(chroms)]
    whole = pd.DataFrame({"chrom": sizes["chrom"], "start": 0, "end": sizes["size"]})
    gpad = genes[["chrom", "start", "end"]].copy()
    gpad["start"] = (gpad["start"] - a.intergenic_gene_pad).clip(lower=0)
    gpad["end"] = gpad["end"] + a.intergenic_gene_pad
    inter_c = subtract(whole, merged_index(gpad, masks))
    for col, val in (("gene_id", "."), ("gene_name", "."), ("strand", "."), ("exon_class", "."),
                     ("max_tpm", np.nan)):
        inter_c[col] = val

    exon_c["compartment"] = "exon"
    intron_c["compartment"] = "intron"
    inter_c["compartment"] = "intergenic"
    comp = pd.concat([exon_c, intron_c, inter_c], ignore_index=True)
    comp = comp[comp["end"] > comp["start"]].sort_values(["chrom", "start"]).reset_index(drop=True)

    # sanity: compartments must not overlap
    for chrom, g in comp.groupby("chrom"):
        s, e = g["start"].to_numpy(), g["end"].to_numpy()
        if np.any(s[1:] < e[:-1]):
            n = int(np.sum(s[1:] < e[:-1]))
            log(f"  WARNING: {n} overlapping compartment intervals on {chrom}; they will be dropped "
                "for TE classification (copies must sit inside a single interval)")

    cols = ["chrom", "start", "end", "compartment", "gene_id", "gene_name", "strand", "exon_class", "max_tpm"]
    comp[cols].to_csv(os.path.join(a.out_dir, "compartments.tsv.gz"), sep="\t", index=False)
    splice_kept.to_csv(os.path.join(a.out_dir, "splice_sites.tsv.gz"), sep="\t", index=False)

    # --- reference transcript per kept gene (Ensembl_canonical if tagged, else longest)
    tx_k = tx[tx["gene_id"].isin(kept["gene_id"])].copy()
    tx_len = kept_ex.assign(L=kept_ex["end"] - kept_ex["start"]).groupby("transcript_id")["L"].sum()
    tx_k["L"] = tx_k["transcript_id"].map(tx_len)
    tx_k = tx_k.sort_values(["gene_id", "canonical", "L"], ascending=[True, False, False])
    ref = tx_k.drop_duplicates("gene_id")
    rex = kept_ex[kept_ex["transcript_id"].isin(ref["transcript_id"])].copy()
    rex["L"] = rex["end"] - rex["start"]
    rex = rex.sort_values(["transcript_id", "start"])
    rex["order"] = np.where(rex["strand"] == "+", rex["start"], -rex["start"])
    rex = rex.sort_values(["transcript_id", "order"])
    rex["tx_offset"] = rex.groupby("transcript_id")["L"].cumsum() - rex["L"]
    rex["tx_len"] = rex["transcript_id"].map(tx_len)
    rex[["gene_id", "transcript_id", "chrom", "strand", "start", "end", "tx_offset", "tx_len"]].to_csv(
        os.path.join(a.out_dir, "ref_tx.tsv.gz"), sep="\t", index=False)

    summ = comp.assign(bp=comp["end"] - comp["start"]).groupby("compartment")["bp"].agg(["count", "sum"])
    log("compartment summary (intervals, bp):\n" + summ.to_string())


if __name__ == "__main__":
    main()
