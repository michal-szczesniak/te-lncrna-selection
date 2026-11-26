#!/usr/bin/env python3
"""Patch 9, step 44. Candidate loci for an illustrative figure: lncRNA splice acceptors provided by antisense Alu
copies (the U-rich complement of the Alu A-rich linker or tail followed by AG).

Input: the per-site table of step 25 (splice_age_sites.tsv.gz), GTEx junction usage of step 26, step-3 copies and
base-to-consensus map (original subfamily consensus), the GENCODE GTF, hg38 and the chain files/genomes of the
other species (same --species arguments as step 25, closest first).

Selection: set lnc_TE_site, Alu, antisense copy, acceptor, preformed (the AG is in the consensus), the junction
seen in GTEx. Two classes, each ranked by GTEx junction reads:
  ape_insertion   the AG is present in the closest species, the copy in the two closest species, and the copy is
                  absent in the two monkeys (4th and 5th species): the exon boundary came with an insertion after
                  the split from Old World monkeys
  lost_in_monkey  the AG is present in the closest species; the copy is present in macaque (4th species) but the AG
                  is absent there
Nothing here enters any statistic of the paper; the output only documents the example shown.

Per candidate (top N per class): gene and transcript models (GTF exons), the Alu copy, GTEx reads of the
junctions of the gene that were analysed, the human sequence around the site in transcript orientation, the
consensus base at each position (original subfamily, transcript orientation), and, per species, the orthologous
base at each position (lifted base by base; '-' if not aligned).

Outputs: <out>/candidates.tsv (one row per candidate), <out>/bundle.json (everything needed to draw the figure)
"""
import argparse
import json
import os

import numpy as np
import pandas as pd
import pysam

from ridl_utils import log, read_gtf, revcomp

COMP = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N", "-": "-"}


class Genome:
    def __init__(self, path):
        if path.endswith(".2bit"):
            import py2bit
            self.tb, self.fa = py2bit.open(path), None
            self.chroms = set(self.tb.chroms().keys())
        else:
            self.fa, self.tb = pysam.FastaFile(path), None
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


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sites", required=True, help="step 25 <out>_sites.tsv.gz")
    ap.add_argument("--usage", required=True, help="step 26 junction_usage.tsv.gz")
    ap.add_argument("--copies", required=True, help="step-3b copies_grouped.tsv.gz")
    ap.add_argument("--map-align", required=True, help="step-3 map_align.parquet (original consensus positions)")
    ap.add_argument("--consensus", required=True, help="step-3 consensus_rec.fa")
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--genome-fasta", required=True)
    ap.add_argument("--species", action="append", default=[], help="name=chain=genome, closest first")
    ap.add_argument("--flank", type=int, default=40, help="nt shown on each side of the acceptor AG")
    ap.add_argument("--top", type=int, default=12)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    S = pd.read_csv(a.sites, sep="\t", low_memory=False)
    names = [c[len("site_"):] for c in S.columns
             if c.startswith("site_") and c not in ("site_id", "site_type") and f"copy_{c[len('site_'):]}" in S.columns]
    log(f"species in the site table: {names}")
    sel = (S["set"] == "lnc_TE_site") & (S["rep_family"] == "Alu") & (S["te_orientation"] == "antisense") \
        & (S["site_type"] == "acceptor") & (S["ancestral"] == "preformed")
    S = S[sel].copy()
    log(f"preformed acceptors in antisense Alu copies: {len(S)}")
    U = pd.read_csv(a.usage, sep="\t")
    U["intron_key"] = U["chrom"] + ":" + U["start"].astype(str) + "-" + U["end"].astype(str) + U["strand"]
    S["intron_key"] = S["site_id"].str[:-1]
    S = S.merge(U[["intron_key", "total_reads", "samples_ge5", "acceptor_share", "max_tissue", "max_tissue_mean"]],
                on="intron_key", how="left")
    S = S[S["total_reads"].fillna(0) > 0]
    log(f"  with GTEx junction reads: {len(S)}")

    ape, mon = names[:2], names[3:5]
    cls = pd.Series("", index=S.index)
    if len(names) >= 5:
        in_ape = S[f"site_{names[0]}"] == "present"
        ins = in_ape & np.logical_and.reduce([S[f"copy_{n}"] == "present" for n in ape]
                                             + [S[f"copy_{n}"] == "absent" for n in mon])
        lost = in_ape & (S[f"copy_{names[3]}"] == "present") & (S[f"site_{names[3]}"] == "absent")
        cls[ins] = "ape_insertion"
        cls[lost & ~ins] = "lost_in_monkey"
    S["locus_class"] = cls
    S = S[S["locus_class"] != ""].sort_values(["locus_class", "total_reads", "acceptor_share"],
                                             ascending=[True, False, False])
    C = S.groupby("locus_class").head(a.top).copy()
    log("candidates: " + ", ".join(f"{k} {v}" for k, v in C["locus_class"].value_counts().items()))
    if C.empty:
        C.to_csv(os.path.join(a.out, "candidates.tsv"), sep="\t", index=False)
        json.dump([], open(os.path.join(a.out, "bundle.json"), "w"))
        return

    chroms = sorted(set(C["chrom"]))
    genes, tx, ex = read_gtf(a.gtf, chroms)
    copies = pd.read_csv(a.copies, sep="\t", dtype={"chrom": str}).set_index("copy_id")
    mp = pd.read_parquet(a.map_align, columns=["copy_id", "pos", "cons_pos"])
    mp = mp[mp["copy_id"].isin(set(C["copy_id"]))]
    cons = load_consensus(a.consensus)
    hg = pysam.FastaFile(a.genome_fasta)

    # positions of the window in transcript order (5'->3' of the RNA); the AG ends at the acceptor
    def window(r):
        lo_, hi_ = int(min(r.pos1, r.pos2)), int(max(r.pos1, r.pos2))
        if r.tx_strand == "+":
            return list(range(lo_ - a.flank, hi_ + 1 + a.flank))
        return list(range(hi_ + a.flank, lo_ - 1 - a.flank, -1))

    def tx_base(b, strand):
        return b if strand == "+" else COMP.get(b, "N")

    wins = {r.site_id: window(r) for r in C.itertuples()}
    lifted = {n: {} for n in names}
    from pyliftover import LiftOver
    for item in a.species:
        name, chain, gpath = item.split("=")
        log(f"{name}: lifting {sum(len(w) for w in wins.values())} bases")
        lo, gen = LiftOver(chain), Genome(gpath)
        for r in C.itertuples():
            seq = []
            for p in wins[r.site_id]:
                m = lo.convert_coordinate(r.chrom, int(p))
                if not m:
                    seq.append("-")
                    continue
                b = gen.base(m[0][0], int(m[0][1]))
                flip = (m[0][2] == "-") != (r.tx_strand == "-")
                seq.append(COMP.get(b, "N") if flip else b)
            lifted[name][r.site_id] = "".join(seq)
        del lo

    bundle, rows = [], []
    for r in C.itertuples():
        w = wins[r.site_id]
        hs = "".join(tx_base(hg.fetch(r.chrom, p, p + 1).upper(), r.tx_strand) for p in w)
        m = mp[mp["copy_id"] == r.copy_id].set_index("pos")["cons_pos"]
        cseq = cons.get(str(r.subfamily_orig), "")
        cp = copies.loc[r.copy_id]
        cb, cpos = [], []
        for p in w:
            q = m.get(p, np.nan)
            if np.isfinite(q) and 0 < q <= len(cseq):
                b = cseq[int(q) - 1]
                cb.append(b if cp["te_strand"] == r.tx_strand else COMP.get(b, "N"))
                cpos.append(int(q))
            else:
                cb.append(" ")
                cpos.append(None)
        g = genes[genes["gene_id"] == r.gene_id]
        gname = g["gene_name"].iloc[0] if len(g) else ""
        t = tx[tx["gene_id"] == r.gene_id]
        e = ex[ex["gene_id"] == r.gene_id]
        models = []
        for tid, te in e.groupby("transcript_id"):
            canon = bool(t.loc[t["transcript_id"] == tid, "canonical"].any()) if len(t) else False
            models.append({"transcript_id": tid, "canonical": canon,
                           "exons": te[["start", "end"]].sort_values("start").values.tolist()})
        gs, ge = (int(g["start"].iloc[0]), int(g["end"].iloc[0])) if len(g) else (int(e["start"].min()), int(e["end"].max()))
        uj = U[(U["chrom"] == r.chrom) & (U["start"] >= gs) & (U["end"] <= ge) & (U["strand"] == r.tx_strand)]
        site = {c: (v.item() if hasattr(v, "item") else v) for c, v in r._asdict().items() if c != "Index"}
        site = {k: (None if isinstance(v, float) and not np.isfinite(v) else v) for k, v in site.items()}
        bundle.append({
            "site": site, "gene_name": gname, "gene_span": [gs, ge], "transcripts": models,
            "copy": {"copy_id": r.copy_id, "chrom": cp["chrom"], "start": int(cp["start"]), "end": int(cp["end"]),
                     "te_strand": cp["te_strand"], "subfamily_orig": str(cp.get("subfamily_orig", "")),
                     "cons_start": int(cp["cons_start"]), "cons_end": int(cp["cons_end"]),
                     "div_pct": float(cp["div_pct"])},
            "junctions": uj[["start", "end", "total_reads", "samples_ge5", "donor_share", "acceptor_share"]]
            .to_dict(orient="records"),
            "window": {"genomic": w, "human": hs, "consensus": "".join(cb), "cons_pos": cpos,
                       "species": {n: lifted[n].get(r.site_id, "") for n in names}},
            "species_order": names,
        })
        rows.append({"locus_class": r.locus_class,
                     "site_id": r.site_id, "gene_id": r.gene_id, "gene_name": gname, "subfamily_orig": r.subfamily_orig,
                     "cons_pos": r.cons_pos, "total_reads": r.total_reads, "samples_ge5": r.samples_ge5,
                     "acceptor_share": r.acceptor_share, "max_tissue": r.max_tissue, "age_class": r.age_class,
                     **{f"copy_{n}": site.get(f"copy_{n}") for n in names},
                     **{f"site_{n}": site.get(f"site_{n}") for n in names},
                     "human_window": hs, "consensus_window": "".join(cb)})
    pd.DataFrame(rows).to_csv(os.path.join(a.out, "candidates.tsv"), sep="\t", index=False)
    with open(os.path.join(a.out, "bundle.json"), "w") as fh:
        json.dump(bundle, fh, indent=1, default=str)
    log(f"wrote {len(rows)} candidates to {a.out}")


if __name__ == "__main__":
    main()
