#!/usr/bin/env python3
"""Step 3b. Group related subfamilies on one coordinate system (e.g. all AluS* on AluSx).

Each member subfamily's consensus (from step 3 --consensus-out) is globally aligned to the group's
reference consensus; consensus positions of every mapped base are translated to reference positions
(bases at member positions that align to a gap in the reference are dropped). Copies of members get
subfamily = group name and keep the original name in subfamily_orig -- use it as a stratum in step 6
(--pop-strata tri,cpg,subfamily_orig --cont-strata subfamily_orig) so that a different subfamily mix
in exons and introns cannot masquerade as a positional signal.

  --group 'AluS*:^AluS:AluSx'      name : regex on subfamily names : reference subfamily
Subfamilies not in any group are passed through unchanged.
Outputs: --out-bases (parquet), --out-copies (tsv.gz), --out-consensus (FASTA incl. group names),
         and a per-member report (identity to the reference, fraction of positions mapped).
"""
import argparse
import re

import numpy as np
import pandas as pd
import pysam

from ridl_utils import log


def nw_map(a, b, match=2, mismatch=-1, gap=-2):
    """Global alignment (linear gaps). Returns (map_a, identity): map_a[i] = 0-based position in b
    aligned to a[i], or -1 if a[i] is aligned to a gap."""
    n, m = len(a), len(b)
    A = np.frombuffer(a.encode(), dtype="S1")
    B = np.frombuffer(b.encode(), dtype="S1")
    G = -gap
    H = np.zeros((n + 1, m + 1), dtype=np.int64)
    P = np.zeros((n + 1, m + 1), dtype=np.int8)  # 0 diag, 1 up, 2 left
    H[0, :] = -G * np.arange(m + 1)
    H[:, 0] = -G * np.arange(n + 1)
    P[0, 1:] = 2
    P[1:, 0] = 1
    jr = np.arange(m + 1)
    for i in range(1, n + 1):
        s = np.where(B == A[i - 1], match, mismatch)
        diag = H[i - 1, :-1] + s
        up = H[i - 1, 1:] - G
        row = np.empty(m + 1, dtype=np.int64)
        row[0] = H[i, 0]
        row[1:] = np.maximum(diag, up)
        ptr = np.zeros(m + 1, dtype=np.int8)
        ptr[1:] = np.where(diag >= up, 0, 1)
        u = row + G * jr
        cm = np.maximum.accumulate(u)
        left = cm > u
        ptr[left] = 2
        ptr[0] = 1
        H[i] = cm - G * jr
        P[i] = ptr
    i, j = n, m
    amap = np.full(n, -1, dtype=np.int64)
    same = cols = 0
    while i > 0 or j > 0:
        p = P[i, j] if i > 0 and j > 0 else (1 if i > 0 else 2)
        if p == 0:
            amap[i - 1] = j - 1
            same += A[i - 1] == B[j - 1]
            cols += 1
            i, j = i - 1, j - 1
        elif p == 1:
            i -= 1
            cols += 1
        else:
            j -= 1
            cols += 1
    return amap, same / max(cols, 1)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bases", required=True, help="mapping parquet from step 3 (align mode)")
    ap.add_argument("--copies", required=True, help="te_copies.tsv.gz from step 2")
    ap.add_argument("--consensus", required=True, help="FASTA from step 3 --consensus-out")
    ap.add_argument("--group", action="append", required=True, help="NAME:REGEX:REFERENCE (repeatable)")
    ap.add_argument("--out-bases", required=True)
    ap.add_argument("--out-copies", required=True)
    ap.add_argument("--out-consensus", required=True)
    ap.add_argument("--report", default=None, help="TSV with identity/mapped fraction per member")
    a = ap.parse_args()

    cons = {}
    with pysam.FastxFile(a.consensus) as fh:
        for rec in fh:
            cons[rec.name.split("#")[0]] = rec.sequence.upper()
    cp = pd.read_csv(a.copies, sep="\t", dtype={"chrom": str})
    cp["subfamily_orig"] = cp["subfamily"]
    present = sorted(cp["subfamily"].unique())

    maps, report, out_cons = {}, [], {}
    assigned = {}
    for spec in a.group:
        name, rx, refname = spec.split(":")
        if refname not in cons:
            raise SystemExit(f"reference {refname} not in consensus FASTA")
        ref = cons[refname]
        out_cons[name] = ref
        for sf in present:
            if not re.search(rx, sf) or sf in assigned:
                continue
            if sf not in cons:
                log(f"  WARNING: no consensus for {sf}; its copies are dropped from group {name}")
                assigned[sf] = None
                continue
            amap, ident = nw_map(cons[sf], ref) if sf != refname else (np.arange(len(ref)), 1.0)
            lut = np.zeros(len(cons[sf]) + 2, dtype=np.int64)  # 1-based member pos -> 1-based ref pos (0 = none)
            lut[1:len(amap) + 1] = np.where(amap >= 0, amap + 1, 0)
            maps[sf] = lut
            assigned[sf] = name
            report.append({"group": name, "subfamily": sf, "reference": refname, "len": len(cons[sf]),
                           "identity": round(float(ident), 3), "frac_mapped": round(float(np.mean(amap >= 0)), 3),
                           "copies": int((cp["subfamily"] == sf).sum())})
    rep = pd.DataFrame(report)
    log("group members:\n" + (rep.to_string(index=False) if len(rep) else "  (none)"))
    if a.report:
        rep.to_csv(a.report, sep="\t", index=False)

    bases = pd.read_parquet(a.bases)
    sf_of = cp.set_index("copy_id")["subfamily"]
    bsf = bases["copy_id"].map(sf_of)
    keep = np.ones(len(bases), dtype=bool)
    newpos = bases["cons_pos"].to_numpy().astype(np.int64).copy()
    for sf, grp in assigned.items():
        m = (bsf == sf).to_numpy()
        if grp is None:
            keep &= ~m
            continue
        lut = maps[sf]
        cpos = np.clip(newpos[m], 0, len(lut) - 1)
        mapped = lut[cpos]
        newpos[m] = mapped
        keep[np.flatnonzero(m)[mapped == 0]] = False
    log(f"  bases kept after remapping: {keep.sum()} / {len(bases)}")
    out = bases.assign(cons_pos=newpos.astype(np.int32))[keep]
    out.to_parquet(a.out_bases, index=False)

    grp_of = {sf: g for sf, g in assigned.items() if g is not None}
    cp = cp[~cp["subfamily"].isin([sf for sf, g in assigned.items() if g is None])].copy()
    cp["subfamily"] = cp["subfamily"].map(lambda x: grp_of.get(x, x))
    cp.to_csv(a.out_copies, sep="\t", index=False)
    with open(a.out_consensus, "w") as fh:
        for name, seq in out_cons.items():
            fh.write(f">{name}\n{seq}\n")
        for sf in present:
            if sf not in assigned and sf in cons:
                fh.write(f">{sf}\n{cons[sf]}\n")
    log("done")


if __name__ == "__main__":
    main()
