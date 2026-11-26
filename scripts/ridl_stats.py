"""Statistics shared by steps 5-7.

Everything is computed on aggregated arrays so that permutations (unit labels) and cluster
bootstrap (gene resampling) only need np.bincount calls.
"""
import glob
import os

import numpy as np
import pandas as pd


def load_annot(path, columns=None, query=None, subfamilies=None):
    """Load per-base annotation (directory of parquet parts or one file).
    columns: list of needed columns (None = all); columns used in `query` are added automatically."""
    import pyarrow.parquet as pq
    files = sorted(glob.glob(os.path.join(path, "*.parquet"))) if os.path.isdir(path) else [path]
    schema = pq.read_schema(files[0]).names
    if columns is not None:
        cols = [c for c in schema if c in set(columns)]
        if query:
            cols += [c for c in schema if c not in cols and c in query]
    else:
        cols = None
    filters = [("subfamily", "in", list(subfamilies))] if subfamilies else None
    parts = []
    for f in files:
        d = pd.read_parquet(f, columns=cols, filters=filters)
        for c in ("subfamily", "compartment", "orientation", "exon_class", "tri", "te_strand"):
            if c in d.columns:
                d[c] = d[c].astype("category")
        parts.append(d)
    df = pd.concat(parts, ignore_index=True)
    for c in ("subfamily", "compartment", "orientation"):
        if c in df.columns:
            df[c] = df[c].astype(str)
    if query:
        df = df.query(query, engine="python")
    return df


def select_comparison(df, subfamily, orientation, fg, bg):
    """Rows for one (subfamily, fg orientation) comparison; label 1 = fg, 0 = bg.
    If the background compartment has no orientation (intergenic), all its copies are used."""
    d = df[df["subfamily"] == subfamily]
    f = d[(d["compartment"] == fg) & (d["orientation"] == orientation)]
    b = d[d["compartment"] == bg]
    if b["orientation"].ne(".").any():
        b = b[b["orientation"] == orientation]
    out = pd.concat([f.assign(_lab=1), b.assign(_lab=0)], ignore_index=True)
    return out


class Agg:
    """Aggregated (unit, cell) table. cell = bin * n_strata + stratum."""

    def __init__(self, d, value_cols, bin_col="bin", strata_col=None, n_bins=None):
        # a unit = one insertion within one label (fragments of an insertion may sit in
        # different compartments; they are treated as separate units then)
        ukey = d["unit_id"].astype(str).to_numpy() + "|" + d["_lab"].astype(str).to_numpy()
        units, uidx = np.unique(ukey, return_inverse=True)
        lab_by_unit = pd.Series(d["_lab"].to_numpy()).groupby(uidx).first().to_numpy()
        gene = d["gene_id"].astype(str).to_numpy()
        clus_key = np.where(gene == ".", d["unit_id"].astype(str).to_numpy(), gene)
        clus_by_unit = pd.Series(clus_key).groupby(uidx).first().to_numpy()
        self.clusters, self.unit_cluster = np.unique(clus_by_unit, return_inverse=True)
        self.units = units
        self.label = lab_by_unit.astype(np.int8)
        b = d[bin_col].to_numpy().astype(np.int64)
        self.n_bins = int(n_bins if n_bins is not None else b.max() + 1)
        if strata_col is None:
            s = np.zeros(len(d), dtype=np.int64)
            self.n_strata = 1
        else:
            s = d[strata_col].to_numpy().astype(np.int64)
            self.n_strata = int(s.max() + 1) if len(s) else 1
        cell = b * self.n_strata + s
        key = uidx.astype(np.int64) * (self.n_bins * self.n_strata) + cell
        uk, inv = np.unique(key, return_inverse=True)
        self.row_unit = (uk // (self.n_bins * self.n_strata)).astype(np.int64)
        self.row_cell = (uk % (self.n_bins * self.n_strata)).astype(np.int64)
        self.vals = {c: np.bincount(inv, weights=d[c].to_numpy().astype(float), minlength=len(uk))
                     for c in value_cols}
        self.vals["_n"] = np.bincount(inv, minlength=len(uk)).astype(float)
        self.n_cells = self.n_bins * self.n_strata

    def sums(self, label, weight):
        """Return dict: name -> (fg array[n_cells], bg array[n_cells])."""
        L = label[self.row_unit]
        w = weight[self.row_unit]
        fgw = w * (L == 1)
        bgw = w * (L == 0)
        out = {}
        for k, v in self.vals.items():
            out[k] = (np.bincount(self.row_cell, weights=v * fgw, minlength=self.n_cells),
                      np.bincount(self.row_cell, weights=v * bgw, minlength=self.n_cells))
        return out

    def units_per_bin(self, label):
        """Distinct units per bin for fg and bg (unweighted)."""
        rb = self.row_cell // self.n_strata
        key = np.unique(self.row_unit * self.n_bins + rb)
        u, b = key // self.n_bins, key % self.n_bins
        fg = np.bincount(b[label[u] == 1], minlength=self.n_bins)
        bg = np.bincount(b[label[u] == 0], minlength=self.n_bins)
        return fg, bg


def mh_logor(a, n1, c, n0, n_bins, n_strata):
    """Mantel-Haenszel log OR per bin and overall. Inputs are per-cell arrays (bin*strata)."""
    b = n1 - a
    d = n0 - c
    N = n1 + n0
    with np.errstate(divide="ignore", invalid="ignore"):
        R = np.where(N > 0, a * d / N, 0.0).reshape(n_bins, n_strata)
        S = np.where(N > 0, b * c / N, 0.0).reshape(n_bins, n_strata)
        per_bin = np.log(R.sum(1) / S.sum(1))
        overall = np.log(R.sum() / S.sum())
    per_bin[~np.isfinite(per_bin)] = np.nan
    return per_bin, (overall if np.isfinite(overall) else np.nan)


def strat_mean_diff(s1, n1, s0, n0, n_bins, n_strata):
    """Difference of means fg-bg per bin, combined across strata with weights n1*n0/(n1+n0)."""
    with np.errstate(divide="ignore", invalid="ignore"):
        m1 = s1 / n1
        m0 = s0 / n0
        h = np.where((n1 > 0) & (n0 > 0), n1 * n0 / (n1 + n0), 0.0)
        diff = np.where(h > 0, m1 - m0, 0.0)
        H = h.reshape(n_bins, n_strata)
        D = (diff * h).reshape(n_bins, n_strata)
        per_bin = D.sum(1) / H.sum(1)
        overall = D.sum() / H.sum()
    per_bin[~np.isfinite(per_bin)] = np.nan
    return per_bin, (overall if np.isfinite(overall) else np.nan)


def permuted_labels(agg, rng, scheme):
    lab = agg.label.copy()
    if scheme == "global":
        return rng.permutation(lab)
    # within cluster (gene): shuffle labels among units of the same gene
    order = np.argsort(agg.unit_cluster, kind="stable")
    cl = agg.unit_cluster[order]
    bounds = np.flatnonzero(np.diff(cl)) + 1
    for grp in np.split(order, bounds):
        if len(grp) > 1:
            lab[grp] = rng.permutation(lab[grp])
    return lab


def bootstrap_weights(agg, rng):
    k = len(agg.clusters)
    draw = np.bincount(rng.integers(0, k, size=k), minlength=k)
    return draw[agg.unit_cluster].astype(float)


def perm_pvalue(obs, null):
    """Two-sided empirical p: null has shape (n_perm, ...)."""
    obs = np.asarray(obs, dtype=float)
    null = np.asarray(null, dtype=float)
    ge = np.sum(np.abs(null) >= np.abs(obs)[None, ...], axis=0)
    valid = np.sum(np.isfinite(null), axis=0)
    p = (1 + ge) / (1 + valid)
    p = np.where(np.isfinite(obs) & (valid > 0), p, np.nan)
    return p


def call_peaks(df, value_col, q_col, q_thr=0.1, min_run=2):
    """Runs of consecutive bins with q < q_thr and the same sign of the effect."""
    peaks = []
    df = df.sort_values("bin")
    sig = (df[q_col] < q_thr).to_numpy() & np.isfinite(df[value_col].to_numpy())
    sign = np.sign(df[value_col].to_numpy())
    bins = df["bin"].to_numpy()
    i = 0
    while i < len(df):
        if not sig[i]:
            i += 1
            continue
        j = i
        while (j + 1 < len(df) and sig[j + 1] and sign[j + 1] == sign[i] and bins[j + 1] == bins[j] + 1):
            j += 1
        if j - i + 1 >= min_run:
            seg = df.iloc[i:j + 1]
            peaks.append({"bin_from": int(bins[i]), "bin_to": int(bins[j]),
                          "cons_from": int(seg["cons_from"].min()), "cons_to": int(seg["cons_to"].max()),
                          "n_bins": j - i + 1, "mean_effect": float(seg[value_col].mean()),
                          "min_q": float(seg[q_col].min())})
        i = j + 1
    return peaks
