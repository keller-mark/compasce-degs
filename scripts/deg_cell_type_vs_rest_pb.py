"""
deg_cell_type_vs_rest_pb.py

Differential expression: one cell type vs. all other cells (pseudobulked, pydeseq2).

Input: pseudobulk H5AD produced by pseudobulk.py with sample_col="all_samples".
  - obs rows are specimens/donors
  - X contains summed raw counts for the target cell type

The comparison is structured as:
  - foreground: specimens from the target cell type (all rows in the input)
  - background: we recompute a "rest" pseudobulk on-the-fly from the full
    normalized H5AD.

Because pydeseq2 expects integer raw counts we work from the 'counts' layer
(summed raw UMI counts), not the log1p-normalized values.

Output CSV columns
------------------
gene            gene symbol / ID
baseMean        mean normalized expression across all specimens
log2FoldChange  log2 fold change (target / rest)
lfcSE           standard error of the LFC
stat            Wald statistic
pvalue          unadjusted p-value
padj            Benjamini-Hochberg adjusted p-value
n_specimens_fg  number of specimens in foreground
n_specimens_bg  number of specimens in background
n_cells_fg_mean mean cells per specimen (foreground)
"""

import argparse
import os
import warnings

import numpy as np
import pandas as pd
import scipy.sparse as sp
from anndata import read_h5ad, AnnData
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

MIN_DONORS = 3
MIN_CELLS_PER_DONOR = 25
DONOR_COL = "patient"
SPECIMEN_COL = "specimen"


def unsanitize(sanitized: str, candidates) -> str:
    for c in candidates:
        if c.replace(" ", "_").replace("/", "_") == sanitized:
            return c
    return sanitized


def load_pseudobulk(pb_path: str) -> AnnData:
    pb = read_h5ad(pb_path)
    # pydeseq2 needs integer counts in X
    if sp.issparse(pb.X):
        X = pb.X.toarray()
    else:
        X = np.array(pb.X)
    pb.X = X.astype(int)
    return pb


def check_feasibility(n_fg: int, n_bg: int) -> bool:
    if n_fg < MIN_DONORS:
        print(f"  SKIP: only {n_fg} foreground specimens (need {MIN_DONORS}).")
        return False
    if n_bg < MIN_DONORS:
        print(f"  SKIP: only {n_bg} background specimens (need {MIN_DONORS}).")
        return False
    return True


def run_deseq2(
    counts_df: pd.DataFrame,
    design_df: pd.DataFrame,
    design_factor: str,
    ref_level: str,
) -> pd.DataFrame:
    """Run pydeseq2 and return a tidy results DataFrame."""
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
        description="Cell-type-vs-rest DEG (pseudobulked, pydeseq2)."
    )
    parser.add_argument("--input", required=True, help="Pseudobulk H5AD (target cell type).")
    parser.add_argument("--output", required=True, help="Output CSV path.")
    parser.add_argument("--cell-type-col", required=True)
    parser.add_argument("--cell-type", required=True)
    # Path to the full normalized H5AD to build the "rest" pseudobulk
    parser.add_argument(
        "--normalized-h5ad",
        required=False,
        default=None,
        help=(
            "Full normalized H5AD path.  Required to compute the 'rest' pseudobulk "
            "on-the-fly.  If omitted the script expects the input pseudobulk H5AD "
            "to already contain both foreground and background obs (group column '_fg')."
        ),
    )
    args = parser.parse_args()

    os.makedirs(os.path.dirname(args.output), exist_ok=True)

    print(f"Loading pseudobulk (foreground): {args.input}")
    fg_pb = load_pseudobulk(args.input)
    n_fg = fg_pb.n_obs
    genes = fg_pb.var_names.tolist()

    fg_counts = pd.DataFrame(
        fg_pb.X, index=fg_pb.obs_names, columns=genes
    )
    fg_counts["_fg"] = 1
    # Record mean cells per specimen for reporting
    n_cells_fg_mean = fg_pb.obs["n_cells"].mean() if "n_cells" in fg_pb.obs.columns else np.nan

    # ------------------------------------------------------------------
    # Build "rest" pseudobulk from full normalized H5AD
    # ------------------------------------------------------------------
    if not args.normalized_h5ad:
        raise ValueError("--normalized-h5ad is required for this script.")

    print(f"Loading full normalized H5AD for 'rest' pseudobulk: {args.normalized_h5ad}")
    full_adata = read_h5ad(args.normalized_h5ad)

    actual_cell_type = unsanitize(
        args.cell_type, full_adata.obs[args.cell_type_col].unique().tolist()
    )
    if actual_cell_type != args.cell_type:
        print(f"  Resolved cell type '{args.cell_type}' -> '{actual_cell_type}'")

    rest_mask = full_adata.obs[args.cell_type_col] != actual_cell_type
    rest_adata = full_adata[rest_mask, :].copy()
    print(f"  'rest' cells: {rest_adata.n_obs}")

    # Sum raw counts per specimen for the "rest" cells
    if "counts" in rest_adata.layers:
        rest_X = rest_adata.layers["counts"]
    else:
        rest_X = rest_adata.X

    rest_obs = rest_adata.obs.copy()
    specimens = rest_obs[SPECIMEN_COL].unique()
    rest_rows = []
    rest_spec_ids = []
    for spec in specimens:
        spec_mask = rest_obs[SPECIMEN_COL] == spec
        n_cells = spec_mask.sum()
        if n_cells < MIN_CELLS_PER_DONOR:
            continue
        spec_X = rest_X[spec_mask.values, :]
        if sp.issparse(spec_X):
            row_sum = np.asarray(spec_X.sum(axis=0)).ravel()
        else:
            row_sum = spec_X.sum(axis=0)
        rest_rows.append(row_sum)
        rest_spec_ids.append(spec)

    if len(rest_rows) == 0:
        raise ValueError("No 'rest' specimens passed the minimum cell threshold.")

    rest_counts = pd.DataFrame(
        np.vstack(rest_rows).astype(int), index=rest_spec_ids, columns=genes
    )

    n_bg = rest_counts.shape[0]
    print(f"  Foreground specimens: {n_fg}, Background specimens: {n_bg}")

    if not check_feasibility(n_fg, n_bg):
        # Write empty output so Snakemake considers the rule done
        pd.DataFrame(
            columns=["gene", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj",
                     "n_specimens_fg", "n_specimens_bg", "n_cells_fg_mean"]
        ).to_csv(args.output, index=False)
        return

    rest_counts["_fg"] = 0

    all_counts = pd.concat([fg_counts, rest_counts])
    design_df = all_counts[["_fg"]].copy()
    design_df["_fg"] = design_df["_fg"].astype(str)
    count_matrix = all_counts.drop(columns=["_fg"])

    print("Running pydeseq2 ...")
    results = run_deseq2(count_matrix, design_df, "_fg", ref_level="0")

    results["n_specimens_fg"] = n_fg
    results["n_specimens_bg"] = n_bg
    results["n_cells_fg_mean"] = round(n_cells_fg_mean, 1) if not np.isnan(n_cells_fg_mean) else np.nan
    results["cell_type"] = args.cell_type
    results["cell_type_col"] = args.cell_type_col

    results.to_csv(args.output, index=False)
    print(f"Wrote {len(results)} rows to {args.output}")


if __name__ == "__main__":
    main()
