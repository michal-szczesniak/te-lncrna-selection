#!/usr/bin/env python3
"""Revision step 38 (patch 5). TEs at lncRNA features compared with sequence right next to the features.

Step 17 (Fig. 4a) compares the share of transcription starts (TSS), splice donors and acceptors and transcript ends
(TES) lying in TEs with the TE share of the exonic sequence of the same genes. Sequence next to a feature can differ
from exonic sequence on average (first exons, exon ends), so here every feature also gets local controls: positions
--dist-min..--dist-max nt away on the exonic side and on the other side of the feature (intron for splice sites,
upstream for TSS, downstream for TES), in the transcript's orientation. Controls are classified by where they fall
(exon / intron of the same gene / outside the gene). Same gene sets, features, RepeatMasker classes and annotation
as step 17 (functions imported from 17_te_architecture.py).

Output: <out>_local.tsv: gene set x feature x control side x distance bin: n and the share in each TE class, with
the share at the features themselves and the ratio feature / control.
"""
import argparse
import importlib.util
import os

import numpy as np
import pandas as pd

from ridl_utils import log, read_gtf


def load_17():
    here = os.path.dirname(os.path.abspath(__file__))
    spec = importlib.util.spec_from_file_location("arch17", os.path.join(here, "17_te_architecture.py"))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def where(pts, ex_merged, genes):
    """'exon' / 'intron' / 'outside' for positions (chrom, pos, gene_id) relative to that gene."""
    out = np.full(len(pts), "outside", dtype=object)
    gb = genes.set_index("gene_id")
    for gid, idx in pts.groupby("gene_id").indices.items():
        p = pts["pos"].to_numpy()[idx]
        s, e = int(gb.at[gid, "start"]), int(gb.at[gid, "end"])
        inside = (p >= s) & (p < e)
        cls = np.where(inside, "intron", "outside").astype(object)
        iv = ex_merged.get(gid)
        if iv is not None and len(iv):
            k = np.searchsorted(iv[:, 0], p, side="right") - 1
            ok = k >= 0
            ok[ok] &= iv[k[ok], 1] > p[ok]
            cls[ok] = "exon"
        out[idx] = cls
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--gene-set", action="append", required=True, help="name=lnc_genes.tsv.gz (step 1)")
    ap.add_argument("--rmsk-out", required=True)
    ap.add_argument("--dist-bins", default="20,50,100,200,500", help="edges of the distance bins (nt)")
    ap.add_argument("--per-feature", type=int, default=2, help="controls per feature, side and distance bin")
    ap.add_argument("--chroms", default=",".join(f"chr{i}" for i in range(1, 23)))
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    m17 = load_17()
    rng = np.random.default_rng(a.seed)
    chroms = a.chroms.split(",")
    genes, tx, ex = read_gtf(a.gtf, chroms)
    rm = m17.read_rmsk(a.rmsk_out, set(chroms))
    edges = [int(x) for x in a.dist_bins.split(",")]
    rows = []
    for item in a.gene_set:
        name, path = item.split("=", 1)
        g = pd.read_csv(path, sep="\t")
        kept = set(g.loc[g["drop_reason"].fillna("") == "", "gene_id"])
        t, e = tx[tx["gene_id"].isin(kept)], ex[ex["gene_id"].isin(kept)]
        f = m17.features(t, e).reset_index(drop=True)
        em = {}
        for gid, d in e.groupby("gene_id"):
            iv = d[["start", "end"]].to_numpy()
            iv = iv[np.argsort(iv[:, 0])]
            mm = [list(iv[0])]
            for s_, e_ in iv[1:]:
                if s_ <= mm[-1][1]:
                    mm[-1][1] = max(mm[-1][1], e_)
                else:
                    mm.append([s_, e_])
            em[gid] = np.array(mm, dtype=np.int64)
        # direction of the exonic side in transcript orientation: TSS and acceptor -> downstream, TES and donor -> upstream
        down = f["feature"].isin(["TSS", "acceptor"]).to_numpy()
        sign_tx = np.where(f["strand"].to_numpy() == "+", 1, -1)
        ctrl = []
        for lo, hi in zip(edges[:-1], edges[1:]):
            for side, sgn in (("exon side", 1), ("other side", -1)):
                for _ in range(a.per_feature):
                    dist = rng.integers(lo, hi, len(f))
                    direction = np.where(down, 1, -1) * sgn * sign_tx
                    c = f[["gene_id", "chrom", "strand", "feature"]].copy()
                    c["pos"] = f["pos"].to_numpy() + direction * dist
                    c["side"], c["dist_bin"] = side, f"{lo}-{hi}"
                    ctrl.append(c)
        C = pd.concat(ctrl, ignore_index=True)
        C["falls_in"] = where(C, em, genes[genes["gene_id"].isin(kept)])
        F = f.assign(side="feature", dist_bin="0", falls_in="feature")
        pts = m17.annotate(pd.concat([F, C], ignore_index=True), rm)
        log(f"{name}: features {len(F)}, controls {len(C)}")
        feat_share = {}
        for (feat, side, db, fi), s in pts.groupby(["feature", "side", "dist_bin", "falls_in"]):
            r = {"gene_set": name, "feature": feat, "side": side, "dist_bin": db, "falls_in": fi, "n": len(s)}
            for cl, fn in m17.CLASSES.items():
                r[cl] = float(np.mean([fn(c, fm) for c, fm in zip(s["rclass"], s["rfamily"])]))
            if side == "feature":
                feat_share[feat] = r
            rows.append(r)
        for r in rows:
            if r["gene_set"] == name and r["side"] != "feature" and r["feature"] in feat_share:
                for cl in m17.CLASSES:
                    v = r[cl]
                    r[f"{cl}_feature_over_control"] = feat_share[r["feature"]][cl] / v if v > 0 else np.nan
    out = pd.DataFrame(rows)
    out.to_csv(f"{a.out}_local.tsv", sep="\t", index=False)
    pd.set_option("display.width", 240)
    show = out[(out["falls_in"].isin(["feature", "exon", "intron"]))][["gene_set", "feature", "side", "dist_bin", "falls_in", "n", "any TE"]]
    log(show.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
