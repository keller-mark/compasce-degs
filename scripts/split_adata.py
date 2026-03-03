from anndata import read_h5ad
import numpy as np
import pandas as pd
import argparse


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-h5ad", type=str, required=True, help = "Path to input H5AD file.")
    parser.add_argument("--output-h5ad", type=str, required=True, help = "Path to output H5AD file.")
    parser.add_argument("--cell-type-col", type=str, required=True, help = "Name of cell type column")
    parser.add_argument("--cell-type-name", type=str, required=True, help = "Cell type to subset for")
    parser.add_argument("--sample-id-col", type=str, required=True, help = "Name of sample ID column")
    parser.add_argument("--sample-id", type=str, required=True, help = "Sample ID to subset for")
    
    args = parser.parse_args()

    adata = read_h5ad(args.input_h5ad)

    cell_type_col = args.cell_type_col
    cell_type_name = args.cell_type_name
    sample_id_col = args.sample_id_col
    sample_id = args.sample_id

    cell_type_mask = adata.obs[cell_type_col] == cell_type_name
    sample_id_mask = adata.obs[sample_id_col] == sample_id
    
    # Subset
    adata_subset = adata[cell_type_mask & sample_id_mask, :].copy()
    
    # Save
    adata_subset.write_h5ad(args.output_h5ad)