"""
deg_within_cell_type_pb.py

Differential expression: within a single cell type, comparing two sample groups
(e.g. AKI vs Healthy Reference), using pseudobulked counts and pydeseq2.

Input: pseudobulk H5AD produced by pseudobulk.py for the relevant sample_col.
  - obs rows are specimens/donors, each annotated with `sample_col`
  - X contains summed raw counts for the target cell type

The wildcard values lhs_group / rhs_group use underscores in place of spaces
(path-sanitized by Snakemake).  We un-sanitize them by matching against the
obs column values, trying both the sanitized and original forms.

Output CSV columns
------------------
gene            gene symbol / ID
baseMean        mean normalized expression across all specimens
log2FoldChange  LHS vs RHS log2 fold change  (positive = up in LHS)
lfcSE           standard error of the LFC
stat            Wald statistic
pvalue          unadjusted p-value
padj            Benjamini-Hochberg adjusted p-value
n_specimens_lhs number of specimens in LHS group
n_specimens_rhs number of specimens in RHS group
n_cells_lhs_mean mean cells per specimen in LHS
n_cells_rhs_mean mean cells per specimen in RHS
"""

import argparse
import os
import warnings

import numpy as np
import pandas as pd
import scipy.sparse as sp
from anndata import read_h5ad
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

MIN_DONORS = 3


def unsanitize(sanitized: str, candidates) -> str:
    """Find the original value from `candidates` that sanitizes to `sanitized`."""
    for c in candidates:
        if c.replace(" ", "_").replace("/", "_") == sanitized:
            return c
    # If not found, return as-is (may already be the original form)
    return sanitized


def check_feasibility(n_lhs: int, n_rhs: int) -> bool:
    if n_lhs < MIN_DONORS:
        print(f"  SKIP: only {n_lhs} LHS specimens (need {MIN_DONORS}).")
        return False
    if n_rhs < MIN_DONORS:
        print(f"  SKIP: only {n_rhs} RHS specimens (need {MIN_DONORS}).")
        return False
    return True


def run_deseq2(
    counts_df: pd.DataFrame,
    design_df: pd.DataFrame,
    design_factor: str,
    ref_level: str,
) -> pd.DataFrame:
    dds = DeseqDataSet(
        counts=counts_df,
        metadata=design_df,
        design_factors=design_factor,
        ref_level=[design_factor, ref_level],
        quiet=True,
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        dds.deseq2()

    stat_res = DeseqStats(dds, quiet=True)
    stat_res.summary()
    results = stat_res.results_df.reset_index().rename(columns={"index": "gene"})
    return results


def main():
    parser = argparse.ArgumentParser(
        description="Within-cell-type DEG between sample groups (pseudobulked, pydeseq2)."
    )
    parser.add_argument("--input", required=True, help="Pseudobulk H5AD path.")
    parser.add_argument("--output", required=True, help="Output CSV path.")
    parser.add_argument("--sample-col", required=True, help="Column name for sample groups (sanitized).")
    parser.add_argument("--lhs-group", required=True, help="LHS group label (sanitized).")
    parser.add_argument("--rhs-group", required=True, help="RHS group label (sanitized).")
    parser.add_argument("--cell-type-col", required=True)
    parser.add_argument("--cell-type", required=True)
    args = parser.parse_args()

    os.makedirs(os.path.dirname(args.output), exist_ok=True)

    print(f"Loading pseudobulk: {args.input}")
    pb = read_h5ad(args.input)

    # Convert counts to dense int for pydeseq2
    if sp.issparse(pb.X):
        X = pb.X.toarray()
    else:
        X = np.array(pb.X)
    X = X.astype(int)

    sample_col_sanitized = args.sample_col
    lhs_sanitized = args.lhs_group
    rhs_sanitized = args.rhs_group

    # Try to find the actual (unsanitized) column name in obs
    actual_col = unsanitize(sample_col_sanitized, pb.obs.columns.tolist())
    if actual_col not in pb.obs.columns:
        raise ValueError(
            f"Column '{actual_col}' (sanitized: '{sample_col_sanitized}') "
            f"not found in pseudobulk obs. Available: {pb.obs.columns.tolist()}"
        )

    # Find actual group labels
    unique_vals = pb.obs[actual_col].unique().tolist()
    actual_lhs = unsanitize(lhs_sanitized, unique_vals)
    actual_rhs = unsanitize(rhs_sanitized, unique_vals)

    print(f"  Comparing '{actual_lhs}' vs '{actual_rhs}' in column '{actual_col}'")

    lhs_mask = pb.obs[actual_col] == actual_lhs
    rhs_mask = pb.obs[actual_col] == actual_rhs

    n_lhs = lhs_mask.sum()
    n_rhs = rhs_mask.sum()
    print(f"  LHS specimens: {n_lhs}, RHS specimens: {n_rhs}")

    n_cells_lhs_mean = pb.obs.loc[lhs_mask, "n_cells"].mean() if "n_cells" in pb.obs.columns else np.nan
    n_cells_rhs_mean = pb.obs.loc[rhs_mask, "n_cells"].mean() if "n_cells" in pb.obs.columns else np.nan

    if not check_feasibility(n_lhs, n_rhs):
        pd.DataFrame(
            columns=["gene", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj",
                     "n_specimens_lhs", "n_specimens_rhs", "n_cells_lhs_mean", "n_cells_rhs_mean"]
        ).to_csv(args.output, index=False)
        return

    # Build count matrix and design for the two groups only
    both_mask = lhs_mask | rhs_mask
    X_sub = X[both_mask.values, :]
    obs_sub = pb.obs[both_mask].copy()
    genes = pb.var_names.tolist()

    counts_df = pd.DataFrame(X_sub, index=obs_sub.index, columns=genes)
    design_df = obs_sub[[actual_col]].copy()
    design_df = design_df.rename(columns={actual_col: "_group"})

    print("Running pydeseq2 ...")
    # ref_level = RHS so that positive LFC means up-regulated in LHS
    results = run_deseq2(counts_df, design_df, "_group", ref_level=actual_rhs)

    results["n_specimens_lhs"] = int(n_lhs)
    results["n_specimens_rhs"] = int(n_rhs)
    results["n_cells_lhs_mean"] = round(n_cells_lhs_mean, 1) if not np.isnan(n_cells_lhs_mean) else np.nan
    results["n_cells_rhs_mean"] = round(n_cells_rhs_mean, 1) if not np.isnan(n_cells_rhs_mean) else np.nan
    results["lhs_group"] = actual_lhs
    results["rhs_group"] = actual_rhs
    results["sample_col"] = actual_col
    results["cell_type"] = args.cell_type
    results["cell_type_col"] = args.cell_type_col

    results.to_csv(args.output, index=False)
    print(f"Wrote {len(results)} rows to {args.output}")


if __name__ == "__main__":
    main()
