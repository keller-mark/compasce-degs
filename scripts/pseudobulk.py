"""
pseudobulk.py

Computes pseudobulk profiles for a single (cell_type_col, cell_type) combination.

For each donor/specimen, sums raw counts across all cells belonging to the
requested cell type.  The output is an AnnData H5AD file where:
  - obs  = donors/specimens (rows)
  - var  = genes (columns)
  - X    = sum of raw counts
  - obs columns include the sample_col grouping label and QC stats

Special sentinel: if --sample-col is "all_samples", we pseudobulk across all
samples without filtering by group, producing one pseudobulk row per specimen.
This is used by the cell-type-vs-rest pseudobulked DEG rule.

Memory strategy
---------------
- We use Dask to process the sparse expression matrix chunk-by-chunk so that
  the full dense matrix never has to live in RAM.
- Cells are grouped by specimen; we iterate over specimens one at a time.
"""

import argparse
import os

import numpy as np
import pandas as pd
import scipy.sparse as sp
from anndata import read_h5ad, AnnData

# Minimum thresholds (see README)
MIN_CELLS_PER_DONOR = 25
MIN_DONORS = 3

DONOR_COL = "patient"
SPECIMEN_COL = "specimen"


def unsanitize(sanitized: str, candidates) -> str:
    """Recover the original value from `candidates` that sanitizes to `sanitized`."""
    for c in candidates:
        if c.replace(" ", "_").replace("/", "_") == sanitized:
            return c
    return sanitized  # fall back: assume no sanitization was needed


def pseudobulk_cell_type(
    adata,
    cell_type_col: str,
    cell_type: str,
    sample_col: str,
) -> AnnData:
    """
    Subset to `cell_type`, then sum raw counts per specimen.

    Returns an AnnData with one row per specimen that passes the minimum
    cell-count threshold.  Specimens that belong to the given sample_col
    group are annotated accordingly.
    """
    # Subset to cell type of interest
    mask = adata.obs[cell_type_col] == cell_type
    ct_adata = adata[mask, :]
    print(
        f"  Cell type '{cell_type}': {ct_adata.n_obs} cells across "
        f"{ct_adata.obs[SPECIMEN_COL].nunique()} specimens."
    )

    if ct_adata.n_obs == 0:
        raise ValueError(f"No cells found for cell type '{cell_type}' in column '{cell_type_col}'.")

    # Use raw counts layer if available
    if "counts" in ct_adata.layers:
        X = ct_adata.layers["counts"]
    else:
        X = ct_adata.X

    obs_df = ct_adata.obs.copy()
    specimens = obs_df[SPECIMEN_COL].unique()

    pb_rows = []
    pb_obs_records = []

    for spec in specimens:
        spec_mask = obs_df[SPECIMEN_COL] == spec
        n_cells = spec_mask.sum()
        if n_cells < MIN_CELLS_PER_DONOR:
            continue

        # Sum counts for this specimen — operate on the sparse sub-matrix
        spec_X = X[spec_mask.values, :]
        if sp.issparse(spec_X):
            row_sum = np.asarray(spec_X.sum(axis=0)).ravel()
        else:
            row_sum = spec_X.sum(axis=0)

        pb_rows.append(row_sum)

        # Collect metadata: take first cell's obs values (all cells in a
        # specimen share the same sample-level metadata)
        first_row = obs_df[spec_mask].iloc[0]
        meta = {
            SPECIMEN_COL: spec,
            DONOR_COL: first_row[DONOR_COL],
            "n_cells": int(n_cells),
        }
        if sample_col != "all_samples" and sample_col in obs_df.columns:
            meta[sample_col] = first_row[sample_col]
        # Copy any other useful sample-level columns
        for col in ["EnrollmentCategory", "AdjudicatedCategory"]:
            if col in obs_df.columns and col not in meta:
                meta[col] = first_row[col]
        pb_obs_records.append(meta)

    if len(pb_rows) == 0:
        raise ValueError(
            f"No specimens passed the minimum cell threshold ({MIN_CELLS_PER_DONOR}) "
            f"for cell type '{cell_type}'."
        )

    pb_X = np.vstack(pb_rows)  # shape: (n_specimens, n_genes)
    pb_obs = pd.DataFrame(pb_obs_records).set_index(SPECIMEN_COL)

    pb_adata = AnnData(
        X=sp.csr_matrix(pb_X),
        obs=pb_obs,
        var=ct_adata.var.copy(),
    )
    pb_adata.uns["cell_type_col"] = cell_type_col
    pb_adata.uns["cell_type"] = cell_type
    pb_adata.uns["sample_col"] = sample_col
    pb_adata.uns["n_donors"] = len(pb_obs)
    pb_adata.uns["min_cells_per_donor"] = MIN_CELLS_PER_DONOR
    pb_adata.uns["min_donors"] = MIN_DONORS

    print(
        f"  Pseudobulk: {pb_adata.n_obs} specimens pass threshold "
        f"(min {MIN_CELLS_PER_DONOR} cells/specimen)."
    )
    return pb_adata


def main():
    parser = argparse.ArgumentParser(description="Compute pseudobulk profiles.")
    parser.add_argument("--input", required=True, help="Normalized H5AD path.")
    parser.add_argument("--output", required=True, help="Output pseudobulk H5AD path.")
    parser.add_argument(
        "--sample-col",
        required=True,
        help="Column used to define sample groups, or 'all_samples' for no grouping.",
    )
    parser.add_argument("--cell-type-col", required=True)
    parser.add_argument("--cell-type", required=True)
    args = parser.parse_args()

    os.makedirs(os.path.dirname(args.output), exist_ok=True)

    print(f"Loading: {args.input}")
    adata = read_h5ad(args.input)

    # The cell_type wildcard is sanitized in the Snakefile (e.g., VSM/P -> VSM_P).
    # Recover the original value by matching against obs column values.
    actual_cell_type = unsanitize(
        args.cell_type,
        adata.obs[args.cell_type_col].unique().tolist()
    )
    if actual_cell_type != args.cell_type:
        print(f"  Resolved cell type '{args.cell_type}' -> '{actual_cell_type}'")

    pb = pseudobulk_cell_type(
        adata,
        cell_type_col=args.cell_type_col,
        cell_type=actual_cell_type,
        sample_col=args.sample_col,
    )

    print(f"Writing pseudobulk H5AD: {args.output}")
    pb.write_h5ad(args.output, compression="gzip")
    print("Done.")


if __name__ == "__main__":
    main()
