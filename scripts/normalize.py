"""
normalize.py

Loads the raw KPMP H5AD file, joins clinical metadata, applies basic
normalization (library-size normalize + log1p) on the counts layer,
then writes the result to a new H5AD file.

Memory notes
------------
- The expression matrix may be >20 GB in dense form; we keep it sparse throughout.
- We write with backed=False (in-memory) since h5ad writing requires the full
  object anyway; the normalized file itself is still sparse CSR/CSC.
"""

import argparse
import os

import numpy as np
import pandas as pd
import scipy.sparse as sp
import scanpy as sc
from anndata import read_h5ad


# ---------------------------------------------------------------------------
# AKI sub-category roll-up: ATI + AIN -> "AKI"
# Per README: "aggregate AKI subcategories ATI and AIN into a parent AKI category"
# ---------------------------------------------------------------------------
AKI_SUBCATEGORIES = {"Acute Tubular Injury", "Acute Interstitial Nephritis"}


def clean_adjudicated_category(row):
    val = row["Primary Adjudicated Category"]
    if pd.notna(val) and val != "NA" and val != "":
        return val
    # Healthy Reference samples are never adjudicated
    if row.get("Enrollment Category") in ("Healthy Reference",):
        return "Healthy Reference"
    return ""


def build_aki_category(row):
    """Roll ATI and AIN up into a parent 'AKI' category."""
    val = row["AdjudicatedCategory"]
    if val in AKI_SUBCATEGORIES:
        return "AKI"
    return val


def load_and_merge(h5ad_path: str, csv_path: str) -> sc.AnnData:
    print(f"Reading H5AD: {h5ad_path}")
    adata = read_h5ad(h5ad_path)

    print(f"Reading clinical CSV: {csv_path}")
    clinical = pd.read_csv(csv_path)

    # Merge on participant ID
    adata.obs = adata.obs.merge(
        clinical, left_on="patient", right_on="Participant ID", how="left"
    )

    # Filter to cells that have clinical data (inner join semantics)
    has_clinical = ~adata.obs["Participant ID"].isna()
    n_before = adata.n_obs
    adata = adata[has_clinical, :].copy()
    print(f"Kept {adata.n_obs}/{n_before} cells after joining clinical data.")

    # Fill NA in Primary Adjudicated Category before applying the clean function
    adata.obs["Primary Adjudicated Category"] = (
        adata.obs["Primary Adjudicated Category"].fillna("NA")
    )

    # Build clean categorical columns
    adata.obs["AdjudicatedCategory"] = adata.obs.apply(
        clean_adjudicated_category, axis="columns"
    )
    adata.obs["EnrollmentCategory"] = adata.obs["Enrollment Category"]

    # Roll AKI sub-categories up
    adata.obs["AdjudicatedCategory"] = adata.obs.apply(
        build_aki_category, axis="columns"
    )

    # Rename dot-separated cell-type columns
    adata.obs = adata.obs.rename(
        columns={
            "subclass.l1": "subclass_l1",
            "subclass.l2": "subclass_l2",
            "subclass.l3": "subclass_l3",
        }
    )

    # Fill remaining NAs in string columns
    for col in adata.obs.columns:
        if pd.api.types.is_string_dtype(adata.obs[col]) or str(adata.obs[col].dtype) == "object":
            adata.obs[col] = adata.obs[col].fillna("NA")

    # Column names cannot contain slashes (Zarr / file-path incompatibility)
    adata.obs = adata.obs.rename(
        columns={c: c.replace("/", " per ") for c in adata.obs.columns}
    )

    return adata


def normalize(adata: sc.AnnData) -> sc.AnnData:
    """Library-size normalise counts and log1p transform.

    We normalise the 'counts' layer (raw UMI counts) if it exists,
    otherwise fall back to adata.X.  The result is stored in adata.X
    and a 'log1p_norm' layer; raw counts are preserved in 'counts'.
    """
    # Identify the raw counts
    if "counts" in adata.layers:
        raw = adata.layers["counts"]
    else:
        raw = adata.X
        adata.layers["counts"] = raw.copy() if sp.issparse(raw) else raw.copy()

    print("Normalising (library-size + log1p) …")
    # Work on a copy so we don't alter the stored raw layer
    adata.X = raw.copy() if sp.issparse(raw) else raw.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.layers["log1p_norm"] = adata.X.copy()

    return adata


def main():
    parser = argparse.ArgumentParser(
        description="Normalize KPMP H5AD and join clinical metadata."
    )
    parser.add_argument("--input-h5ad", required=True)
    parser.add_argument("--input-csv", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    os.makedirs(os.path.dirname(args.output), exist_ok=True)

    adata = load_and_merge(args.input_h5ad, args.input_csv)
    adata = normalize(adata)

    print(f"Writing normalized H5AD to: {args.output}")
    adata.write_h5ad(args.output, compression="gzip")
    print("Done.")


if __name__ == "__main__":
    main()
