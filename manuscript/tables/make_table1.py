#!/usr/bin/env python3
"""Table 1 of the manuscript ("What the data exclude, and what the test can see") from the pipeline outputs.

Columns and sources:
  genes in the gene-level test      23_bounds.py, table main, both orientations, layer excess_poly_noncpg (genes)
  largest TE-specific reduction     the same row, bound_pct (lower 95% bound relative to the intronic frequency)
  largest TE-specific gain          upper 95% bound of excess_phylop_p447 / excess_phylop_prim
  maximal share of exonic copies    bound_pct / (0.2 x 0.5): a domain over a fifth of the element losing half of
                                    its variants (Methods: phi = limit / (d r))
  conserved TE stretches            24_te_positive_control.py, all compartments, comparison A: Alu groups with the
                                    top 1% of windows of primate phyloP, older groups with mammalian phyloP above 1.0
  planted 20-nt domain              33_spikein_scan.py with genes resampled: power_window for half of the variants
                                    removed in half of the exonic copies, median over positions, sense / antisense
Usage: python tables/make_table1.py --project <dir> [--rev ...] [--rev5b ...] [--check manuscript.md] --out <md>
"""
import argparse
from decimal import Decimal, ROUND_HALF_UP
import os
import re

import numpy as np
import pandas as pd

GROUPS = [("AluS*", "AluS", "alu_prim_q0.99", "AluS"), ("AluJ*", "AluJ", "alu_prim_q0.99", "AluJ"),
          ("MIR*", "MIR", "old_p447_t1.0", "MIR"), ("L2abc", "L2a/b/c", "old_p447_t1.0", "L2abc"),
          ("LTR", "LTR", None, None), ("L1", "L1", None, None)]
HEAD = ["TE group", "Genes in the gene-level test",
        "Largest TE-specific reduction in variant frequency consistent with the data (95%)",
        "Largest TE-specific gain in conservation consistent with the data (phyloP, mammals / primates, 95%)",
        "Maximal share of exonic copies with a domain over 20% of the element losing half of its variants",
        "Conserved TE stretches: fewer human variants",
        "Planted 20-nt domain, half of variants, half of copies: detection, median over positions (sense / antisense)"]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--project", required=True)
    ap.add_argument("--rev", default=None)
    ap.add_argument("--rev5b", default=None)
    ap.add_argument("--check", default=None, help="manuscript markdown with Table 1")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    rev = a.rev or os.path.join(a.project, "results", "revision")
    rev5b = a.rev5b or os.path.join(a.project, "results", "revision5b")
    B = pd.read_csv(os.path.join(rev, "bounds_bounds.tsv"), sep="\t")
    B = B[(B.table == "main") & (B.orientation == "pooled")]
    LB = {"LTR": pd.read_csv(os.path.join(rev, "ltr_bounds_class_bounds.tsv"), sep="\t"),
          "L1": pd.read_csv(os.path.join(rev, "l1lin_bounds_class_bounds.tsv"), sep="\t")}
    rows = []
    for g, lab, pc, scan in GROUPS:
        t = LB[g][LB[g].orientation == "pooled"] if g in LB else B[B.subfamily == g]
        t = t.set_index("layer")
        v = t.loc["excess_poly_noncpg"]
        gain = f"{t.loc['excess_phylop_p447', 'ci_hi']:.2f} / {t.loc['excess_phylop_prim', 'ci_hi']:.2f}"
        share = f"{v.bound_pct / (0.2 * 0.5):.0f}%" if pc else "–"
        cons = "–"
        if pc:
            p = pd.read_csv(os.path.join(rev, f"posctrl_{pc}_posctrl.tsv"), sep="\t")
            x = p[(p.group == g) & (p.compartment == "all") & p.comparison.str.startswith("A")].iloc[0]
            cons = f"{100 * (1 - np.exp(x.mh_log_or)):.1f}%"
        plant = "–"
        if scan:
            med = []
            for o in ("sense", "antisense"):
                s = pd.read_csv(os.path.join(rev5b, f"{scan}_{o}_power.tsv"), sep="\t")
                s = s[(s.segment_nt == 20) & np.isclose(s.reduction, 0.5) & np.isclose(s.fraction_units, 0.5)]
                # medians of detection shares are often x.xx5; round half up, as in the manuscript
                q = Decimal(str(round(float(np.median(s.power_window)), 6)))
                med.append(str(q.quantize(Decimal("0.01"), rounding=ROUND_HALF_UP)))
            plant = " / ".join(med)
        rows.append([lab, f"{int(v.genes)}", f"{v.bound_pct:.1f}%", gain, share, cons, plant])
    md = ["| " + " | ".join(HEAD) + " |", "| " + " | ".join("---" for _ in HEAD) + " |"] + \
        ["| " + " | ".join(r) + " |" for r in rows]
    open(a.out, "w").write("\n".join(md) + "\n")
    print("\n".join(md))
    if a.check:
        ref = open(a.check).read()
        sec = ref[ref.index("## Table 1"):]
        tab = [l for l in sec.split("\n") if l.startswith("|")][2:2 + len(rows)]
        n = 0
        for r, l in zip(rows, tab):
            c = [x.strip() for x in l.strip().strip("|").split("|")]
            for x, y in zip(r, c):
                if x != y:
                    n += 1
                    print(f"  differs: '{x}' vs manuscript '{y}'")
        print(f"cells differing from the manuscript: {n}")


if __name__ == "__main__":
    main()
