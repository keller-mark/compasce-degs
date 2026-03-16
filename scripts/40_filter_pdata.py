# TODO: filter the pseudobulked .pdata file based on the number of cells and number of samples thresholds.

# These thresholds will be provided via argparse params via the snakefile.

# agg_adata.py sets the "num_cells_orig" column which we can use for the number of cells filtering.

from anndata import read_h5ad, AnnData
import numpy as np
import pandas as pd
import argparse

NUM_CELLS_COLNAME = "num_cells_orig"

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-h5ad", type=str, required=True, help = "Path to input H5AD file.")
    parser.add_argument("--output-h5ad", type=str, required=True, help = "Path to output H5AD file.")
    parser.add_argument("--cell-type-col", type=str, required=True, help = "Name of cell type column")
    parser.add_argument("--sample-id-col", type=str, required=True, help = "Name of sample ID column")
    parser.add_argument("--num-samples-threshold", type=int, required=True, help = "Min number of samples when using pseudobulked data")
    parser.add_argument("--num-cells-per-sample-threshold", type=int, required=True, choices=["mean", "sum"], default="sum", help = "Min number of cells per sample when using pseudobulked data")
    
    args = parser.parse_args()

    adata = read_h5ad(args.input_h5ad)