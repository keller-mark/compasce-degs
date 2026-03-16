from anndata import read_h5ad, AnnData
import numpy as np
import pandas as pd
import pertpy as pt
import argparse

NUM_CELLS_COLNAME = "num_cells_orig"

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-h5ad", type=str, required=True, help = "Path to input H5AD file.")
    parser.add_argument("--output-de-csv", type=str, required=True, help = "Path to output CSV file for DE results.")
    parser.add_argument("--output-filtering-csv", type=str, required=True, help = "Path to output CSV file for filtering results, for debugging.")
    parser.add_argument("--cell-type-col", type=str, required=True, help = "Name of cell type column")
    parser.add_argument("--sample-id-col", type=str, required=True, help = "Name of sample ID column")
    parser.add_argument("--cell-type-name", type=str, required=True, help = "Cell type to subset for")
    parser.add_argument("--num-samples-threshold", type=int, required=True, help = "Min number of samples per group threshold")
    parser.add_argument("--num-cells-per-sample-threshold", type=int, required=True, help = "Min number of cells per sample threshold")
    
    args = parser.parse_args()

    pdata = read_h5ad(args.input_h5ad)

    print(f"Orig pdata shape: {pdata.shape}")

    # Filter pdata
    sample_id_col = args.sample_id_col
    cell_type = args.cell_type_name
    cell_type_col = args.cell_type_col
    num_samples_threshold = args.num_samples_threshold
    num_cells_per_sample_threshold = args.num_cells_per_sample_threshold

    # First filter based on number of cells per sample. We can use the "num_cells_orig" column that we set in agg_adata.py to do this filtering.

    # These (cell type, sample ID) rows have at least `num_cells_per_sample_threshold` cells, so we keep them.
    pdata.obs["has_sufficient_num_cells"] = pdata.obs[NUM_CELLS_COLNAME] >= num_cells_per_sample_threshold

    # For each cell type, count the number of unique sample IDs that have sufficient number of cells,
    # and filter out cell types that don't have at least `num_samples_threshold` samples with sufficient number of cells.
    
    num_samples_per_cell_type = (
        pdata.obs[pdata.obs["has_sufficient_num_cells"]]
            .groupby(cell_type_col)[sample_id_col]
            .nunique()
    )

    pdata.obs["has_sufficient_num_samples"] = pdata.obs.apply(lambda row: num_samples_per_cell_type[row[cell_type_col]] >= num_samples_threshold, axis=1)

    # Save this filtering info and write for debugging
    filtering_by_num_cells_df = pdata.obs[["has_sufficient_num_cells", "has_sufficient_num_samples", cell_type_col, sample_id_col, NUM_CELLS_COLNAME]].copy()
    filtering_by_num_cells_df.to_csv(args.output_filtering_csv, index=True)
    
    # Subset the anndata object.
    pdata = pdata[pdata.obs["has_sufficient_num_cells"] & pdata.obs["has_sufficient_num_samples"]].copy()

    print(f"Filtered pdata shape: {pdata.shape}")
    
    # Cell type vs rest
    pdata.obs["is_cell_type"] = pdata.obs[cell_type_col].apply(lambda x: cell_type if x == cell_type else "rest")

    # For cell type vs. rest, we use a design such as "~subclass_l1"
    # Reference: https://hbctraining.github.io/DGE_workshop/lessons/04_DGE_DESeq2_analysis.html
    pds2 = pt.tl.PyDESeq2(adata=pdata, design=f"~is_cell_type")
    pds2.fit()

    df = pds2.test_contrasts(pds2.contrast(column="is_cell_type", baseline="rest", group_to_compare=cell_type))

    print("Done with DE analysis, writing output CSV...")

    df.to_csv(args.output_de_csv, index=True)
