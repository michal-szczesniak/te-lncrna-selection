#!/usr/bin/env python3
"""Revision step 31. Build --group arguments for 03b_group_subfamilies.py from family patterns, choosing the
reference automatically and keeping only close relatives.

For each NAME:REGEX given with --family, the subfamilies present in the step-2 copies that match REGEX are
collected; the member with most copies becomes the reference; each other member's reconstructed consensus
(step 3 --consensus-out) is globally aligned to the reference (the aligner of step 3b) and joins the group only
if its identity is >= --min-identity. Members below the threshold stay separate subfamilies (as MIR3 and minor
L2 families in the paper). Prints the arguments one token per line ('--group', then 'NAME:^(A|B|...)$:REF'),
to be read into a bash array (mapfile -t G < <(python 31_make_groups.py ...); ... "${G[@]}"), and writes a
report with identities to --report. Use group names without shell wildcards.
"""
import argparse
import importlib.util
import os
import re
import sys

import pandas as pd
import pysam


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--copies", required=True)
    ap.add_argument("--consensus", required=True)
    ap.add_argument("--family", action="append", required=True, help="NAME:REGEX")
    ap.add_argument("--min-identity", type=float, default=0.75)
    ap.add_argument("--min-copies", type=int, default=50)
    ap.add_argument("--report", required=True)
    a = ap.parse_args()
    here = os.path.dirname(os.path.abspath(__file__))
    spec = importlib.util.spec_from_file_location("g3b", os.path.join(here, "03b_group_subfamilies.py"))
    g3b = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(g3b)
    cons = {}
    with pysam.FastxFile(a.consensus) as fh:
        for rec in fh:
            cons[rec.name.split("#")[0]] = rec.sequence.upper()
    cp = pd.read_csv(a.copies, sep="\t", usecols=["subfamily"])
    n = cp["subfamily"].value_counts()
    rows, args, used = [], [], set()
    for spec_ in a.family:
        name, rx = spec_.split(":", 1)
        members = [s for s in n.index if re.search(rx, s) and s in cons and s not in used and n[s] >= a.min_copies]
        if not members:
            print(f"  no members for {name}", file=sys.stderr)
            continue
        ref = max(members, key=lambda s: n[s])
        keep = []
        for s in members:
            ident = 1.0 if s == ref else g3b.nw_map(cons[s], cons[ref])[1]
            rows.append({"group": name, "subfamily": s, "reference": ref, "copies": int(n[s]),
                         "identity": round(float(ident), 3), "joined": ident >= a.min_identity})
            if ident >= a.min_identity:
                keep.append(s)
        used.update(keep)
        if len(keep) >= 1:
            args += ["--group", f"{name}:^({'|'.join(re.escape(k) for k in keep)})$:{ref}"]
    pd.DataFrame(rows).to_csv(a.report, sep="\t", index=False)
    print("\n".join(args))


if __name__ == "__main__":
    main()
