"""
deg_cell_type_vs_rest_nonpb.py

Differential expression: one cell type vs. all other cells (not pseudobulked).

Uses scanpy's `rank_genes_groups` with the Wilcoxon rank-sum test, which is the
closest equivalent to Seurat's FindAllMarkers (wilcox method).

We run the test directly on the single-cell level (no pseudobulking), keeping the
expression matrix sparse throughout to avoid OOM errors.

Output CSV columns
------------------
gene            gene name / ID
score           Wilcoxon z-score
log2fc          log2 fold change (log1p_norm scale)
pval            uncorrected p-value
pval_adj        Benjamini-Hochberg adjusted p-value
n_cells_fg      number of cells in the foreground (cell type of interest)
n_cells_bg      number of cells in the background (all other cells)
n_donors_fg     number of unique donors in the foreground
n_donors_bg     number of unique donors in the background
"""

import argparse
import os
import warnings

import numpy as np
import pandas as pd
import scanpy as sc
from anndata import read_h5ad

DONOR_COL = "patient"


def unsanitize(sanitized: str, candidates) -> str:
    for c in candidates:
        if c.replace(" ", "_").replace("/", "_") == sanitized:
            return c
    return sanitized


def run_deg(adata: sc.AnnData, cell_type_col: str, cell_type: str) -> pd.DataFrame:
    # Binary grouping: foreground = cell type of interest, background = rest
    fg_mask = adata.obs[cell_type_col] == cell_type
    bg_mask = ~fg_mask

    n_fg = fg_mask.sum()
    n_bg = bg_mask.sum()
    n_donors_fg = adata.obs.loc[fg_mask, DONOR_COL].nunique() if DONOR_COL in adata.obs.columns else np.nan
    n_donors_bg = adata.obs.loc[bg_mask, DONOR_COL].nunique() if DONOR_COL in adata.obs.columns else np.nan

    print(f"  Foreground: {n_fg} cells ({n_donors_fg} donors), Background: {n_bg} cells ({n_donors_bg} donors)")

    if n_fg == 0:
        raise ValueError(f"No cells found for cell type '{cell_type}' in column '{cell_type_col}'.")

    # Assign binary labels for rank_genes_groups
    adata.obs["_group"] = "rest"
    adata.obs.loc[fg_mask, "_group"] = "target"

    # Use log1p_norm layer if available (normalized + log-transformed), else X
    if "log1p_norm" in adata.layers:
        adata_tmp = adata.copy()
        adata_tmp.X = adata_tmp.layers["log1p_norm"]
    else:
        adata_tmp = adata

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sc.tl.rank_genes_groups(
            adata_tmp,
            groupby="_group",
            groups=["target"],
            reference="rest",
            method="wilcoxon",
            use_raw=False,
            pts=True,
        )

    result = sc.get.rank_genes_groups_df(adata_tmp, group="target")
    result = result.rename(
        columns={
            "names": "gene",
            "scores": "score",
            "logfoldchanges": "log2fc",
            "pvals": "pval",
            "pvals_adj": "pval_adj",
        }
    )

    result["n_cells_fg"] = int(n_fg)
    result["n_cells_bg"] = int(n_bg)
    result["n_donors_fg"] = int(n_donors_fg) if not np.isnan(n_donors_fg) else np.nan
    result["n_donors_bg"] = int(n_donors_bg) if not np.isnan(n_donors_bg) else np.nan
    result["cell_type"] = cell_type
    result["cell_type_col"] = cell_type_col

    # Drop the temporary grouping column from adata.obs
    adata.obs.drop(columns=["_group"], inplace=True, errors="ignore")

    return result


def main():
    parser = argparse.ArgumentParser(
        description="Cell-type-vs-rest DEG (not pseudobulked, Wilcoxon)."
    )
    parser.add_argument("--input", required=True, help="Normalized H5AD path.")
    parser.add_argument("--output", required=True, help="Output CSV path.")
    parser.add_argument("--cell-type-col", required=True)
    parser.add_argument("--cell-type", required=True)
    args = parser.parse_args()

    os.makedirs(os.path.dirname(args.output), exist_ok=True)

    print(f"Loading: {args.input}")
    adata = read_h5ad(args.input)

    actual_cell_type = unsanitize(
        args.cell_type, adata.obs[args.cell_type_col].unique().tolist()
    )
    if actual_cell_type != args.cell_type:
        print(f"  Resolved cell type '{args.cell_type}' -> '{actual_cell_type}'")

    print(f"Running Wilcoxon DEG for '{actual_cell_type}' vs rest ...")
    df = run_deg(adata, args.cell_type_col, actual_cell_type)

    df.to_csv(args.output, index=False)
    print(f"Wrote {len(df)} rows to {args.output}")


if __name__ == "__main__":
    main()
