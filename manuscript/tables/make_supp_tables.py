#!/usr/bin/env python3
"""Supplemental Tables S1-S9 from the pipeline outputs.

Every number in the supplemental tables is computed here from result files of the pipeline; the captions are the
text of the manuscript, kept in supplement_template.md, in which each table is a placeholder ({{TABLE:<id>}}).
Tables S8 (ENCODE accessions) and S9 (data sources) are lists, not results, and are taken from the template.

Usage:
  python tables/make_supp_tables.py --project <project dir> [--rev results/revision] [--rev5 results/revision5]
         [--rev5b results/revision5b] [--rev6 results/revision6] [--template tables/supplement_template.md]
         [--check <reference markdown>] --out <out dir>

Inputs (relative to --project unless given): work/attrition.tsv and work/attrition_genes.tsv (step 20, S1);
work/prof_*/<group>_b20*/ and work/comp*/ (steps 6 and 14, S3B, S3D, S3F); work/splice/ltr_tss.tsv (S5A);
results/revision (patches 1-4), results/revision5 (patch 5), results/revision5b (patch 5b), results/revision6
(patch 6, optional: canonical-exon limits and well-supported gene subsets are added as columns when present).

Outputs: <out>/supplemental_tables.md (all tables with captions), <out>/S<n>.md (one file per table),
<out>/check.tsv (with --check: every regenerated cell compared with the reference, e.g. the manuscript's tables).
"""
import argparse
import os
import re
import sys

import numpy as np
import pandas as pd

LAB = {"AluS*": "AluS", "AluJ*": "AluJ", "AluY*": "AluY", "MIR*": "MIR", "MIR3": "MIR3", "L2abc": "L2a/b/c",
       "L2d*": "L2d", "LTR": "LTR", "L1": "L1", "AluS": "AluS", "AluJ": "AluJ", "AluY": "AluY", "MIR": "MIR",
       "nonTE": "repeat-free"}
MISSING = []


def need(path):
    if not os.path.exists(path):
        MISSING.append(path)
        raise FileNotFoundError(path)
    return path


def rd(path):
    return pd.read_csv(need(path), sep="\t", low_memory=False)


def sgn(x, d=3, sign=True):
    """number with a true minus sign; '–' for missing"""
    if x is None or (isinstance(x, float) and not np.isfinite(x)) or pd.isna(x):
        return "–"
    if abs(x) < 0.5 * 10 ** -d:
        x = 0.0
    s = f"{x:+.{d}f}" if sign else f"{x:.{d}f}"
    return s.replace("-", "−")


def ci(e, lo, hi, d=3):
    return f"{sgn(e, d)} [{sgn(lo, d, False)}, {sgn(hi, d, False)}]"


def table(header, rows):
    return ["| " + " | ".join(header) + " |", "| " + " | ".join("---" for _ in header) + " |"] + \
        ["| " + " | ".join(str(c) for c in r) + " |" for r in rows]


# ======================================================================== S1 selection of genes and exonic sequence
def t_S1(P):
    t = rd(os.path.join(P.work, "attrition.tsv"))
    lab = {"exons of kept lncRNA genes": "Exons of kept genes",
           "minus splice sites": "Margins around splice sites (±10 nt)",
           "minus transcription starts": "Margins around transcription starts (±500 nt)",
           "minus transcript ends": "Margins around transcript ends (±100 nt)",
           "minus segmental duplications": "Segmental duplications",
           "minus other genes' spans": "Spans of other genes"}
    rows = []
    for _, r in t.iterrows():
        if r.step not in lab:
            continue
        rows.append([lab[r.step], "–" if r.removed_alone_bp == 0 else f"{r.removed_alone_bp / 1e6:.2f} Mb",
                     f"{r.left_bp / 1e6:.2f} Mb", "100%" if r.left_share == 1 else f"{100 * r.left_share:.1f}%",
                     f"{int(r.genes_with_exons):,}"])
    return {"S1": table(["Step", "Removed alone", "Left", "Share left", "Genes with exonic sequence"], rows)}


# ======================================================================== S2 LTR and L1 groups
def t_S2(P):
    rows = []
    for cls, f, fm in (("LTR", "groups_auto_ltr.tsv", "groups_ltr.tsv"), ("L1", "groups_auto_l1.tsv", "groups_l1.tsv")):
        t = rd(os.path.join(P.rev, f))
        mp = rd(os.path.join(P.rev, fm)).set_index("subfamily")
        for _, r in t.iterrows():
            cons = r.reference if r.joined else f"own ({r.subfamily})"
            mapped = f"{100 * mp.loc[r.subfamily, 'frac_mapped']:.1f}%" if r.joined and r.subfamily in mp.index else "–"
            rows.append([cls, r.group.replace("_g", ""), r.subfamily, f"{int(r.copies):,}", f"{100 * r.identity:.1f}%",
                         cons, mapped])
    return {"S2": table(["Class", "Family pattern", "Subfamily", "Copies", "Identity to reference", "Consensus used",
                         "Mapped"], rows)}


# ======================================================================== S3A limits
def t_S3A(P):
    B = rd(os.path.join(P.rev, "bounds_bounds.tsv"))
    B = B[B.layer == "excess_poly_noncpg"]
    LB = {"LTR": rd(os.path.join(P.rev, "ltr_bounds_class_bounds.tsv")),
          "L1": rd(os.path.join(P.rev, "l1lin_bounds_class_bounds.tsv"))}
    C6 = None
    f6 = os.path.join(P.rev6, "canonical_bounds_bounds.tsv") if P.rev6 else None
    if f6 and os.path.exists(f6):
        C6 = pd.read_csv(f6, sep="\t")
        C6 = C6[(C6.layer == "excess_poly_noncpg") & (C6.table == "canonical")]

    def bound(g, tbl="main", ori="pooled", w=False):
        if g in LB:
            if tbl != "main":
                return None
            t = LB[g]
            t = t[(t.layer == "excess_poly_noncpg") & (t.orientation == ori)]
        elif tbl == "canonical":
            if C6 is None:
                return None
            t = C6[(C6.subfamily == g) & (C6.orientation == ori)]
        else:
            tb = f"embedded_{g.rstrip('*')}" if tbl == "embedded" else tbl
            t = B[(B.table == tb) & (B.subfamily == g) & (B.orientation == ori)]
        if not len(t):
            return None
        r = t.iloc[0]
        return f"{(r.w_bound_pct if w else r.bound_pct):.1f}% ({int(r.genes)})"

    head = ["TE group", "Main", "Weighted", "Sense only", "Antisense only", "Embedded copies", "Extended gene set",
            "100-nt margin"] + (["Canonical-transcript exons"] if C6 is not None else [])
    rows = []
    for g in ("AluS*", "AluJ*", "MIR*", "L2abc", "LTR", "L1"):
        cells = [bound(g), bound(g, w=True), bound(g, ori="sense"), bound(g, ori="antisense"), bound(g, "embedded"),
                 bound(g, "extended"), bound(g, "tss100")] + ([bound(g, "canonical")] if C6 is not None else [])
        rows.append([LAB[g]] + [c if c else "–" for c in cells])
    return {"S3A": table(head, rows)}


# ======================================================================== S3B whole-element tests in five analyses
def t_S3B(P):
    W = P.work
    runs = [("Main", "prof_v5/{g}_b20_paired", "comp/te_{g}_paired"),
            ("Main, within gene", "prof_gene/{g}_b20", None),
            ("100-nt margin", "prof_tss100/{g}_b20_paired", "comp_tss100/te_{g}_paired"),
            ("Extended", "prof_ext/{g}_b20_paired", "comp_ext/te_{g}_paired"),
            ("Extended, within gene", "prof_ext/{g}_b20_gene", None)]
    rows, n_found = [], 0
    for g in ("AluS", "AluJ", "MIR", "L2abc"):
        cells = {k: [] for k in ("vs", "va", "ps", "pa", "gc", "pk")}
        for tag, prof, comp in runs:
            f = os.path.join(W, prof.format(g=g), "whole_element.tsv")
            if os.path.exists(f):
                n_found += 1
                w = pd.read_csv(f, sep="\t")

                def one(layer, ori):
                    x = w[(w.layer == layer) & (w.orientation == ori)]
                    # values rounded to 4 decimals first, as in the robustness table of the frozen pipeline
                    # (compare_checks.py), then to 3 for display
                    r4 = lambda v: round(float(v), 4)
                    return ci(r4(x.effect.iloc[0]), r4(x.ci_lo.iloc[0]), r4(x.ci_hi.iloc[0])) if len(x) else "–"
                cells["vs"].append(one("population", "sense"))
                cells["va"].append(one("population", "antisense"))
                cells["ps"].append(one("phylop_p447", "sense"))
                cells["pa"].append(one("phylop_p447", "antisense"))
                pkf = os.path.join(W, prof.format(g=g), "peaks_centered.tsv")
                try:
                    cells["pk"].append(str(len(pd.read_csv(pkf, sep="\t"))))
                except (FileNotFoundError, pd.errors.EmptyDataError):
                    cells["pk"].append("0")
            else:
                MISSING.append(f)
                for k in ("vs", "va", "ps", "pa", "pk"):
                    cells[k].append("–")
            s = None
            if comp:
                fs = os.path.join(W, comp.format(g=g) + "_summary.tsv")
                if os.path.exists(fs):
                    t = pd.read_csv(fs, sep="\t").set_index(["metric", "comparison"])
                    s = t.loc[("phylop_p447", "gc_matched")]
            cells["gc"].append(ci(round(s.exon_minus_intron, 4), round(s.ci_lo, 4), round(s.ci_hi, 4))
                               if s is not None else "–")
        for k, lab in (("vs", "Variants, sense (log OR)"), ("va", "Variants, antisense (log OR)"),
                       ("ps", "phyloP, sense"), ("pa", "phyloP, antisense"), ("gc", "phyloP, GC-matched"),
                       ("pk", "Centered peaks")):
            rows.append([LAB[g], lab] + cells[k])
    if not n_found:
        raise FileNotFoundError("no whole-element profiles")
    return {"S3B": table(["TE group", "Measure"] + [r[0] for r in runs], rows)}


# ======================================================================== S3C GC matching on copy or flanks
def t_S3C(P):
    def fg(s, mtr, comp):
        t = rd(os.path.join(P.rev, f"flankgc_{s}_paired_summary.tsv"))
        r = t[(t.metric == mtr) & (t.comparison == comp)].iloc[0]
        return ci(r.exon_minus_intron, r.ci_lo, r.ci_hi)
    rows = []
    for s, n in (("pc", "protein-coding, repeat-free"), ("lnc", "lncRNA, repeat-free"), ("AluS", "lncRNA, AluS"),
                 ("AluJ", "lncRNA, AluJ"), ("MIR", "lncRNA, MIR"), ("L2abc", "lncRNA, L2a/b/c"),
                 ("LTR", "lncRNA, LTR"), ("L1", "lncRNA, L1")):
        for mtr, tr in (("phylop_p447", "mammals"), ("phylop_prim", "primates")):
            rows.append([n, tr] + [fg(s, mtr, c) for c in ("unmatched", "gc_own", "gc_flank", "gc_both")])
    return {"S3C": table(["Set", "phyloP", "Unmatched", "GC of the unit", "GC of the flanks", "Both"], rows)}


# ======================================================================== S3D, S3E, S3F (patch-4 parts)
def _m4(x, d=3):
    if pd.isna(x):
        return "–"
    s = f"{x:+.{d}f}" if d else f"{x:+.0f}"
    if float(s) == 0:
        s = f"{0:.{d}f}"
    return s.replace("-", "−")


def _ci4(e, lo, hi, d=3):
    return f"{_m4(e, d)} [{_m4(lo, d).lstrip('+')}, {_m4(hi, d).lstrip('+')}]"


def t_S3D(P):
    rows = []
    for g in ["AluS", "AluJ", "AluY", "MIR", "MIR3", "L2abc"]:
        d = rd(os.path.join(P.work, "prof_v5", f"{g}_b20", "whole_element.tsv"))
        for o in ["sense", "antisense"]:
            p = d[(d.layer == "population") & (d.orientation == o)].iloc[0]
            s = d[(d.layer == "singletons") & (d.orientation == o)].iloc[0]
            rows.append([LAB[g], o, f"{int(p.units_fg):,}", _ci4(p.effect, p.ci_lo, p.ci_hi),
                         _ci4(s.effect, s.ci_lo, s.ci_hi)])
    nt = os.path.join(P.work, "prof_nonte", "exon_vs_intron", "whole_element.tsv")
    if os.path.exists(nt):  # repeat-free 100-nt windows of the same compartments, for comparison
        d = rd(nt)
        p, s = (d[d.layer == k].iloc[0] for k in ("population", "singletons"))
        rows.append(["Repeat-free", "–", f"{int(p.units_fg):,} windows", _ci4(p.effect, p.ci_lo, p.ci_hi),
                     _ci4(s.effect, s.ci_lo, s.ci_hi)])
    return {"S3D": table(["TE group", "Orientation", "Exonic units", "Variants (log OR)", "Singleton share"], rows)}


def t_S3E(P):
    t = rd(os.path.join(P.rev, "dropout_te_genes.tsv"))
    rows = []
    f = lambda x, y, d=1: f"{x:.{d}f} / {y:.{d}f}".replace("-", "−")
    for g in ["AluS*", "AluJ*", "AluY*", "MIR*", "MIR3", "L2abc"]:
        for o in ["sense", "antisense"]:
            a = t[(t.group == g) & (t.orientation == o) & (t.paired == False)].iloc[0]  # noqa: E712
            b = t[(t.group == g) & (t.orientation == o) & (t.paired == True)].iloc[0]  # noqa: E712
            rows.append([LAB[g], o, f"{int(a.genes)} / {int(b.genes)}", f(a.median_span / 1000, b.median_span / 1000),
                         f"{100 * a.share_single_exon:.0f}% / {100 * b.share_single_exon:.0f}%",
                         f(a.median_splice, b.median_splice, 2), f(a.median_synth, b.median_synth, 2),
                         f(a.median_max_tpm, b.median_max_tpm, 1)])
    return {"S3E": table(["TE group", "Orientation", "Genes (exonic only / both)", "Gene span, kb", "Single-exon genes",
                          "Splicing efficiency", "Synthesis rate", "Maximum TPM"], rows)}


def t_S3F(P):
    rows = []
    for f, lab in [("lnc_paired_bins.tsv", "repeat-free"), ("te_AluS_paired_bins.tsv", "AluS"),
                   ("te_AluJ_paired_bins.tsv", "AluJ"), ("te_MIR_paired_bins.tsv", "MIR"),
                   ("te_L2abc_paired_bins.tsv", "L2a/b/c")]:
        d = rd(os.path.join(P.work, "comp", f))
        e = d.units_exon / d.units_exon.sum()
        i = d.units_intron / d.units_intron.sum()
        rows.append([lab, f"{int(d.units_exon.sum()):,}", f"{int(d.units_intron.sum()):,}",
                     f"{100 * e.iloc[-2:].sum():.0f}%", f"{100 * i.iloc[-2:].sum():.0f}%"])
    return {"S3F": table(["Set", "Exonic units", "Intronic units", "Exonic units in the top two GC quintiles",
                          "Intronic units in the top two GC quintiles"], rows)}


# ======================================================================== S3G limits in gene subsets
def _m1(x, d=1):
    if pd.isna(x):
        return "–"
    s = f"{x:+.{d}f}"
    if float(s) == 0:
        s = f"{0:.{d}f}"
    return s.replace("-", "−")


def _ci1(e, lo, hi, d=1):
    return f"{_m1(e, d)} [{_m1(lo, d).lstrip('+')}, {_m1(hi, d).lstrip('+')}]"


def t_S3G(P):
    b = rd(os.path.join(P.rev5, "bounds_sub_bounds.tsv"))
    x = b[(b.layer == "excess_poly_noncpg") & (b.orientation == "pooled")]
    T = [("all", "All genes"), ("expr_top25", "Top quarter by expression"), ("expr_top10", "Top tenth by expression"),
         ("broad_tau25", "Broadest quarter (tissue specificity τ)")]
    f6 = os.path.join(P.rev6, "bounds_well_bounds.tsv") if P.rev6 else None
    if f6 and os.path.exists(f6):
        b6 = pd.read_csv(f6, sep="\t")
        b6 = b6[(b6.layer == "excess_poly_noncpg") & (b6.orientation == "pooled") & (b6.table != "all")]
        x = pd.concat([x, b6], ignore_index=True)
        T += [(k, l) for k, l in (("gtex_any", "Genes with GTEx estimate"), ("lincrna_v29", "lincRNAs of GENCODE 29"),
                                  ("lincrna_v29_gtex", "GENCODE 29 lincRNAs with GTEx estimate"))
              if (b6.table == k).any()]
        # a subset with the same genes as another one (e.g. every GENCODE 29 lincRNA has a GTEx estimate) is not repeated
        sub6 = os.path.join(P.rev6, "bounds_well_subsets.tsv")
        if os.path.exists(sub6):
            n6 = pd.read_csv(sub6, sep="\t").set_index("subset")["genes"]
            if {"lincrna_v29", "lincrna_v29_gtex"} <= set(n6.index) and n6["lincrna_v29"] == n6["lincrna_v29_gtex"]:
                T = [(k, l) for k, l in T if k != "lincrna_v29_gtex"]
    rows = []
    for g in ["AluS*", "AluJ*", "MIR*", "L2abc"]:
        cells = []
        for k, _ in T:
            r = x[(x.subfamily == g) & (x.table == k)]
            if len(r) == 0 or pd.isna(r.iloc[0].mean_pct):
                cells.append("–")
                continue
            r = r.iloc[0]
            cells.append(f"{_ci1(r.mean_pct, r.ci_lo_pct, r.ci_hi_pct)} ({int(r.genes)}); limit {r.bound_pct:.1f}%")
        rows.append([LAB[g]] + cells)
    return {"S3G": table(["TE group"] + [l for _, l in T], rows)}


# ======================================================================== S3H, S3I comparisons without exclusions
def t_S3H(P):
    def get(v, g, d):
        f = os.path.join(P.rev5, f"{v}_{g}_{d}_summary.tsv")
        if os.path.exists(f):
            return pd.read_csv(f, sep="\t").set_index(["metric", "comparison"]), f
        if v == "main" and d == "paired":
            f = os.path.join(P.rev, f"flankgc_{'lnc' if g == 'nonTE' else g}_paired_summary.tsv")
            if os.path.exists(f):
                return pd.read_csv(f, sep="\t").set_index(["metric", "comparison"]), f
        MISSING.append(f)
        return None, None

    def units(f):
        b = f.replace("_summary", "_bins")
        if not os.path.exists(b) or "flankgc" in f:
            return None
        b = pd.read_csv(b, sep="\t")
        return int(b.units_exon.sum()), int(b.units_intron.sum())
    cols = [("main", "all"), ("tss0", "all"), ("open", "all"), ("main", "paired"), ("tss0", "paired"),
            ("open", "paired")]
    cl = {"main": "main", "tss0": "no start/end margins", "open": "no exclusions"}
    rows = []
    for g in ["AluS", "AluJ", "AluY", "MIR", "MIR3", "L2abc", "nonTE"]:
        cells = []
        for v, d in cols:
            s, f = get(v, g, d)
            u = units(f) if f else None
            cells.append(f"{u[0]:,} / {u[1]:,}" if u else "–")
        rows.append([LAB[g], "units (exon / intron)"] + cells)
        for comp, lab in [("unmatched", "unmatched"), ("gc_both", "GC-matched")]:
            cells = []
            for v, d in cols:
                s, f = get(v, g, d)
                if s is None:
                    cells.append("–")
                    continue
                x = s.loc[("phylop_p447", comp)]
                cells.append(f"{_m4(x.exon_minus_intron)} [{_m4(x.ci_lo).lstrip('+')}, {_m4(x.ci_hi).lstrip('+')}]")
            rows.append([LAB[g], lab] + cells)
    return {"S3H": table(["Group", "Matching"] + [f"{cl[v]}, {'all genes' if d == 'all' else 'same genes'}"
                                                  for v, d in cols], rows)}


def t_S3I(P):
    def f2(x, lo, hi):
        if pd.isna(x) or not np.isfinite(x):
            return "–"
        s = f"{x:.2f}"
        if pd.notna(lo) and pd.notna(hi):
            s += f" [{lo:.2f}, {hi:.2f}]"
        return s
    d = rd(os.path.join(P.rev5b, "conserved_pooled.tsv"))
    d = d[(d.track == "p447") & (d.threshold == 1.5)]
    V = [("main", "main"), ("tss0", "no start/end margins"), ("open", "no exclusions")]
    out = {}
    for des, key in (("all genes", "S3I_all"), ("same genes", "S3I_same")):
        head = ["Group"]
        for _, l in V:
            head += [f"{l}: conserved units (exon / intron)", f"{l}: odds ratio", f"{l}: ratio to repeat-free"]
        rows = []
        for g in ["AluS*", "AluJ*", "AluY*", "MIR*", "MIR3", "L2abc", "nonTE"]:
            cells = []
            for v, _ in V:
                r = d[(d.variant == v) & (d.group == g) & (d.design == des)]
                if len(r) == 0:
                    cells += ["–"] * 3
                    continue
                r = r.iloc[0]
                cells.append(f"{int(r.conserved_exon):,} / {int(r.conserved_intron):,} (of {int(r.units_exon):,} / "
                             f"{int(r.units_intron):,})")
                ok = r.conserved_exon > 0 and r.conserved_intron > 0
                cells.append(f2(r.mh_or, r.ci_lo, r.ci_hi) if ok else "–")
                cells.append(f2(r.ratio_vs_nonte, r.ratio_lo, r.ratio_hi) if (ok and g != "nonTE") else "–")
            rows.append([LAB[g]] + cells)
        out[key] = table(head, rows)
    return out


# ======================================================================== S4 positive control
def t_S4(P):
    pos = [("alu_prim_q0.99", "primates", "top 1% of windows"), ("alu_prim_q0.95", "primates", "top 5% of windows"),
           ("alu_prim_t1.0", "primates", "mean > 1.0"), ("old_p447_t1.0", "mammals", "mean > 1.0"),
           ("old_p447_t1.5", "mammals", "mean > 1.5"), ("old_p447_q0.99", "mammals", "top 1% of windows")]
    rows = []
    pct = lambda v: 100 * (1 - np.exp(v))
    for f, track, thr in pos:
        t = rd(os.path.join(P.rev, f"posctrl_{f}_posctrl.tsv"))
        for g in t.group.unique():
            for comp in ("all", "exon", "intron"):
                x = t[(t.group == g) & (t.compartment == comp)]
                if not len(x):
                    continue
                A = x[x.comparison.str.startswith("A")].iloc[0]
                Bb = x[x.comparison.str.startswith("B")].iloc[0]
                if A.units * A.share_units_with_stretch < 10:
                    continue
                if not np.isfinite(A.mh_log_or):
                    ca = cb = "–"
                else:
                    ca = f"{pct(A.mh_log_or):.1f}% [{pct(A.ci_hi):.1f}, {pct(A.ci_lo):.1f}]".replace("-", "−")
                    cb = f"{pct(Bb.mh_log_or):.1f}% [{pct(Bb.ci_hi):.1f}, {pct(Bb.ci_lo):.1f}]".replace("-", "−")
                sing = (f"{sgn(A.singleton_diff)} [{sgn(A.singleton_ci_lo, 3, False)}, {sgn(A.singleton_ci_hi, 3, False)}]"
                        if np.isfinite(A.singleton_diff) else "–")
                rows.append([LAB.get(g, g), f"{track}, {thr}", comp, f"{int(A.units):,}",
                             f"{100 * A.share_units_with_stretch:.1f}%", ca, cb, sing])
    return {"S4": table(["TE group", "phyloP, threshold", "Compartment", "Units", "Share with stretch",
                         "Fewer variants (A)", "Fewer variants (B)", "Singleton difference (A)"], rows)}


# ======================================================================== S5 TEs at features
def t_S5A(P):
    t = rd(os.path.join(P.work, "splice", "ltr_tss.tsv")).sort_values("lnc_TSS", ascending=False).head(25)
    rows = [[r["name"], r.rfamily, f"{int(r.lnc_TSS_n):,}", f"{100 * r.lnc_TSS:.2f}%", f"{100 * r.lnc_bg_exon:.2f}%",
             f"{r.lnc_TSS_obs_exp:.2f}", f"{100 * r.pc_TSS:.3f}%", f"{r.lnc_vs_pc_TSS:.0f}"] for _, r in t.iterrows()]
    return {"S5A": table(["Subfamily", "Family", "lncRNA starts", "Share of lncRNA starts",
                          "Share of lncRNA exonic sequence", "Observed/expected", "Share of protein-coding starts",
                          "lncRNA/protein-coding"], rows)}


def t_S5B(P):
    d = rd(os.path.join(P.rev5, "feature_local_local.tsv"))
    feat = [("TSS", "transcription starts"), ("donor", "splice donors"), ("acceptor", "splice acceptors"),
            ("TES", "transcript ends")]

    def share(z, k):
        return (z[k] * z.n).sum() / z.n.sum()
    p = lambda v: f"{100 * v:.1f}%"
    rows = []
    for gs, gl in [("lnc", "lncRNA"), ("pc", "protein-coding")]:
        for ft, fl in feat:
            y = d[(d.gene_set == gs) & (d.feature == ft)]
            f = y[y.side == "feature"].iloc[0]
            e1 = y[(y.side == "exon side") & (y.falls_in == "exon") & (y.dist_bin.isin(["20-50", "50-100"]))]
            e2 = y[(y.side == "exon side") & (y.falls_in == "exon") & (y.dist_bin.isin(["100-200", "200-500"]))]
            ot = y[(y.side == "other side") & (y.dist_bin.isin(["20-50", "50-100"]))]
            rows.append([gl, fl, f"{int(f.n):,}", p(f["any TE"]), p(share(e1, "any TE")), p(share(e2, "any TE")),
                         p(share(ot, "any TE")), f"{f['any TE'] / share(e2, 'any TE'):.2f}",
                         f"{p(f['Alu'])} / {p(share(e1, 'Alu'))}", f"{p(f['LTR'])} / {p(share(e1, 'LTR'))}"])
    return {"S5B": table(["Genes", "Feature", "Features", "Any TE: feature", "exonic side 20–100 nt",
                          "exonic side 100–500 nt", "other side 20–100 nt", "ratio (feature / 100–500 nt)", "Alu: feature / exonic side",
                          "LTR: feature / exonic side"], rows)}


def t_S5C(P):
    bn = rd(os.path.join(P.rev5, "boundary_null_summary.tsv"))
    rows = []
    for g in ["AluJ", "AluS", "AluY", "MIR", "MIR3", "L2abc"]:
        for ori in ["antisense", "sense"]:
            r = bn[(bn.group == g) & (bn.orientation == ori)].iloc[0]
            rows.append([LAB[g], ori, f"{int(r.units):,}", f"{int(r.genes):,}", f"{int(r.median_copy_len)}",
                         f"{100 * r.observed:.0f}% [{100 * r.obs_lo:.0f}, {100 * r.obs_hi:.0f}]",
                         f"{100 * r.expected:.0f}% [{100 * r.exp_lo:.0f}, {100 * r.exp_hi:.0f}]", f"{r.ratio:.2f}"])
    return {"S5C": table(["TE group", "Orientation", "Copies", "Genes", "Median copy length (nt)", "Observed",
                          "Expected", "Observed/expected"], rows)}


# ======================================================================== S6 splice-site age
SP = [("panTro6", "chimpanzee"), ("gorGor6", "gorilla"), ("ponAbe3", "orangutan"), ("rheMac10", "macaque"),
      ("calJac4", "marmoset"), ("mm39", "mouse")]


def t_S6A(P):
    A = rd(os.path.join(P.rev, "splice_age_summary.tsv"))
    subs = [("lnc_TE_site", "all", "TE-derived, all"), ("lnc_TE_site", "ancestral=preformed", "TE-derived, preformed"),
            ("lnc_TE_site", "ancestral=created", "TE-derived, created"),
            ("lnc_TE_site", "site_type=acceptor", "TE-derived, acceptors"),
            ("lnc_TE_site", "site_type=donor", "TE-derived, donors"),
            ("lnc_TE_site", "rep_family=Alu", "TE-derived, Alu"), ("lnc_TE_site", "rep_family=MIR", "TE-derived, MIR"),
            ("lnc_TE_site", "rep_family=L2", "TE-derived, L2"),
            ("lnc_TE_site", "te_orientation=antisense", "TE-derived, antisense copies"),
            ("lnc_TE_site", "te_orientation=sense", "TE-derived, sense copies"),
            ("lnc_nonTE_site", "all", "lncRNA, outside repeats"), ("pc_site", "all", "protein-coding")]

    def share(s, sub, meas, d=1):
        r = A[(A.set == s) & (A.subset == sub) & (A.measure == meas)]
        return f"{100 * r.iloc[0].share:.{d}f}" if len(r) else "–"
    rows = []
    for s, sub, lab in subs:
        r = A[(A.set == s) & (A.subset == sub)]
        n = int(r.sites.iloc[0]) if len(r) else 0
        pres = [share(s, sub, f"present_in_{sp} (of aligned)") for sp, _ in SP]
        young = 100 * sum(float(A[(A.set == s) & (A.subset == sub) & (A.measure == f"age={k}")].share.sum())
                          for k in ("human-specific", "shared with panTro6", "shared with gorGor6",
                                    "shared with ponAbe3"))
        pre = share(s, sub, "preformed") if s == "lnc_TE_site" else "–"
        rows.append([lab, f"{n:,}", pre] + pres + [f"{young:.1f}"])
    return {"S6A": table(["Sites", "n", "Preformed (%)"] + [f"Present in {n} (%)" for _, n in SP] +
                         ["After the macaque split (%)"], rows)}


def t_S6B(P):
    b = rd(os.path.join(P.rev, "splice_age_base_baseline.tsv"))
    sets = [("lnc_TE_site", "all", "TE-derived, all"), ("lnc_TE_site", "ancestral=preformed", "TE-derived, preformed"),
            ("lnc_TE_site", "ancestral=created", "TE-derived, created"), ("lnc_TE_site", "rep_family=Alu", "TE-derived, Alu"),
            ("lnc_TE_site", "rep_family=MIR", "TE-derived, MIR"), ("lnc_TE_site", "rep_family=L2", "TE-derived, L2"),
            ("lnc_nonTE_site", "all", "lncRNA, outside repeats"), ("pc_site", "all", "protein-coding")]
    rows = []
    for st, sub, lab in sets:
        for cond in (["all aligned", "copy present"] if st == "lnc_TE_site" else ["all aligned"]):
            cells = []
            for sp, _ in SP:
                x = b[(b.set == st) & (b.subset == sub) & (b.condition == cond) & (b.species == sp)]
                if len(x) == 0:
                    cells.append("–")
                    continue
                x = x.iloc[0]
                cells.append(f"{_m4(100 * x['diff'], 1)} [{_m4(100 * x.diff_lo, 1).lstrip('+')}, "
                             f"{_m4(100 * x.diff_hi, 1).lstrip('+')}] ({100 * x.ctrl_present:.0f}%)")
            rows.append([lab, "copy shared" if cond == "copy present" else "all aligned"] + cells)
    return {"S6B": table(["Sites", "Condition"] + [n for _, n in SP], rows)}


def t_S6C(P):
    O = rd(os.path.join(P.rev, "splice_age_summary.tsv"))
    A = rd(os.path.join(P.rev5b, "splice_age_ca_summary.tsv"))
    ages = ["age=human-specific", "age=shared with panTro6", "age=shared with gorGor6", "age=shared with ponAbe3"]
    subs = [("all", "TE-derived, all"), ("rep_family=Alu", "Alu"), ("rep_family=MIR", "MIR"), ("rep_family=L2", "L2"),
            ("ancestral=preformed", "preformed"), ("ancestral=created", "created")]

    def sh(T, sub, m):
        x = T[(T.set == "lnc_TE_site") & (T.subset == sub)].set_index("measure")
        return x.share.get(m, 0.0), int(x.sites.iloc[0])
    rows = []
    for sub, lab in subs:
        ao = sum(sh(O, sub, m)[0] for m in ages)
        aa = sum(sh(A, sub, m)[0] for m in ages)
        mo, n = sh(O, sub, "age=shared with mm39")
        ma = sh(A, sub, "age=shared with mm39")[0]
        rows.append([lab, f"{n:,}", f"{100 * ao:.1f}%", f"{100 * aa:.1f}%", f"{100 * mo:.1f}%", f"{100 * ma:.1f}%"])
    return {"S6C": table(["Subset", "Sites", "After the macaque split: site rule", "After the macaque split: copy rule",
                          "Present in mouse (most distant): site rule", "Present in mouse: copy rule"], rows)}


# ======================================================================== S7 splice-site constraint by stratum
def t_S7A(P):
    S = rd(os.path.join(P.rev, "splice_strata_total_reads_strata.tsv"))
    cmp_ = [("protein-coding splice sites vs intronic flank", "protein-coding vs intronic bases"),
            ("lncRNA non-TE splice sites vs intronic flank", "lncRNA outside repeats vs intronic bases"),
            ("lncRNA TE-derived splice sites vs unused copies (main)", "TE-derived vs unused copies"),
            ("lncRNA TE-derived splice sites vs intronic", "TE-derived vs intronic copies"),
            ("lncRNA TE-derived splice sites vs unused copies (main) (strata x CpG)",
             "TE-derived vs unused copies, CpG strata"),
            ("lncRNA TE-derived splice sites vs intronic (strata x CpG)", "TE-derived vs intronic copies, CpG strata")]
    sub_ = [("all", "all"), ("no reads", "no GTEx reads"), ("usage q1", "usage, low"), ("usage q2", "usage, middle"),
            ("usage q3", "usage, high"), ("ancestral=preformed", "preformed"), ("ancestral=created", "created"),
            ("ancestral=created&rep_family=Alu", "created, Alu"), ("ancestral=created&rep_family=MIR", "created, MIR"),
            ("ancestral=created&rep_family=L2", "created, L2")]
    # patch 16: pooled rows of Fig. 7A for the two sets outside TEs (main splice-site analysis, 200 resamples),
    # which the usage-strata table does not contain
    M = rd(os.path.join(P.work, "splice", "res_stats.tsv"))
    rows = []
    for c, cl in cmp_:
        if "TE-derived" not in cl and not len(S[(S.comparison == c) & (S.subset == "all")]):
            x = M[(M.comparison == c) & (M.subset == "all")]
            if len(x):
                r = x.iloc[0]
                rows.append([cl, "all (Fig. 7A)", f"{int(r.bases_site):,}",
                             ci(r.mh_log_or, r.mh_log_or_lo, r.mh_log_or_hi, 2),
                             ci(r.singleton_diff, r.singleton_diff_lo, r.singleton_diff_hi, 3),
                             ci(r.phylop_p447_diff, r.phylop_p447_diff_lo, r.phylop_p447_diff_hi, 2)])
        for sub, sl in sub_:
            x = S[(S.comparison == c) & (S.subset == sub)]
            if not len(x) or ("TE-derived" not in cl and sub.startswith("ancestral")):
                continue
            r = x.iloc[0]
            rows.append([cl, sl, f"{int(r.bases_site):,}", ci(r.mh_log_or, r.mh_log_or_lo, r.mh_log_or_hi, 2),
                         ci(r.singleton_diff, r.singleton_diff_lo, r.singleton_diff_hi, 3),
                         ci(r.phylop_p447_diff, r.phylop_p447_diff_lo, r.phylop_p447_diff_hi, 2)])
    return {"S7A": table(["Comparison", "Subset", "Bases", "Variants (log OR)", "Singleton share", "phyloP (mammals)"],
                         rows)}


def t_S7B(P):
    s = rd(os.path.join(P.rev5, "usage_summary_sets.tsv")).set_index("set")
    t = rd(os.path.join(P.rev5, "usage_summary_tertiles.tsv"))
    rows = []
    for k, lab in [("lnc_TE_site", "TE-derived lncRNA"), ("lnc_nonTE_site", "lncRNA outside repeats"),
                   ("pc_site", "protein-coding")]:
        r = s.loc[k]
        z = t[t.set == k].sort_values("tertile")
        f = lambda q: f"{int(q.reads_min):,}–{int(q.reads_max):,} ({int(q.reads_median):,})"
        rows.append([lab, f"{int(r.sites):,}", f"{100 * r.no_reads:.0f}%"] + [f(q) for _, q in z.iterrows()] +
                    [" / ".join(f"{int(q.samples_ge1_median):,}" for _, q in z.iterrows())])
    return {"S7B": table(["Sites", "Sites analyzed", "Without usage estimate", "Low: reads", "Middle: reads",
                          "High: reads", "Samples with reads (low / middle / high)"], rows)}


ORDER = ["S1", "S2", "S3A", "S3B", "S3C", "S3D", "S3E", "S3F", "S3G", "S3H", "S3I_all", "S3I_same", "S4", "S5A",
         "S5B", "S5C", "S6A", "S6B", "S6C", "S7A", "S7B", "S8", "S9"]
BUILDERS = [t_S1, t_S2, t_S3A, t_S3B, t_S3C, t_S3D, t_S3E, t_S3F, t_S3G, t_S3H, t_S3I, t_S4, t_S5A, t_S5B, t_S5C,
            t_S6A, t_S6B, t_S6C, t_S7A, t_S7B]


def cells(line):
    return [c.strip().replace("\\[", "[").replace("\\]", "]") for c in line.strip().strip("|").split("|")]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--project", required=True)
    ap.add_argument("--rev", default=None)
    ap.add_argument("--rev5", default=None)
    ap.add_argument("--rev5b", default=None)
    ap.add_argument("--rev6", default=None)
    ap.add_argument("--work", default=None, help="default <project>/work")
    ap.add_argument("--template", default=os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                                       "supplement_template.md"))
    ap.add_argument("--check", default=None, help="reference markdown with the same tables, e.g. the manuscript's")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    P = argparse.Namespace(work=a.work or os.path.join(a.project, "work"),
                           rev=a.rev or os.path.join(a.project, "results", "revision"),
                           rev5=a.rev5 or os.path.join(a.project, "results", "revision5"),
                           rev5b=a.rev5b or os.path.join(a.project, "results", "revision5b"),
                           rev6=a.rev6 or os.path.join(a.project, "results", "revision6"))
    os.makedirs(a.out, exist_ok=True)
    tpl = open(a.template).read()
    gen = {}
    for b in BUILDERS:
        try:
            gen.update(b(P))
        except (FileNotFoundError, KeyError, IndexError) as e:
            print(f"[skip] {b.__name__}: {type(e).__name__} {e}", file=sys.stderr)
    # fill the template: generated tables replace placeholders; tables not generated keep the template's
    # stored copy (<!-- keep:ID --> blocks), which is reported
    out = tpl
    kept = []
    for k in ORDER:
        ph = "{{TABLE:" + k + "}}"
        if ph not in out:
            continue
        stored = re.search(r"<!-- keep:" + k + r" -->\n(.*?)\n<!-- /keep -->", out, flags=re.S)
        if k in gen:
            body = "\n".join(gen[k])
        elif stored:
            body = stored.group(1)
            if k not in ("S8", "S9"):
                kept.append(k)
        else:
            body = f"[table {k} not generated: inputs missing]"
            kept.append(k)
        out = out.replace(ph, body)
    out = re.sub(r"<!-- keep:\w+ -->\n.*?\n<!-- /keep -->\n?", "", out, flags=re.S)
    open(os.path.join(a.out, "supplemental_tables.md"), "w").write(out)
    for part in re.split(r"\n(?=### Supplement(?:al|ary) Table )", out):
        m = re.match(r"### Supplement(?:al|ary) Table (S\d)", part)
        if m:
            open(os.path.join(a.out, f"{m.group(1)}.md"), "w").write(part.strip() + "\n")
    print(f"generated: {', '.join(k for k in ORDER if k in gen)}")
    if kept:
        print(f"NOT regenerated (template copy kept): {', '.join(kept)}")
    if MISSING:
        print("missing inputs:\n  " + "\n  ".join(sorted(set(MISSING))))
    if a.check:
        ref = open(a.check).read()
        refs = re.findall(r"((?:^\|.*\n?)+)", ref, flags=re.M)
        res = []
        for k in ORDER:
            if k not in gen or ORDER.index(k) >= len(refs):
                continue
            l1 = [cells(x) for x in gen[k][2:]]
            l2 = [cells(x) for x in refs[ORDER.index(k)].strip().split("\n")[2:]]
            n_diff = sum(c1 != c2 for r1, r2 in zip(l1, l2) for c1, c2 in zip(r1, r2))
            n_diff += abs(len(l1) - len(l2)) + sum(abs(len(r1) - len(r2)) for r1, r2 in zip(l1, l2))
            res.append({"table": k, "rows": len(l1), "rows_ref": len(l2), "cells_differing": n_diff})
            for r1, r2 in zip(l1, l2):
                for c1, c2 in zip(r1, r2):
                    if c1 != c2:
                        print(f"  {k}: '{c1}' vs reference '{c2}'")
        pd.DataFrame(res).to_csv(os.path.join(a.out, "check.tsv"), sep="\t", index=False)
        print(pd.DataFrame(res).to_string(index=False))


if __name__ == "__main__":
    main()
