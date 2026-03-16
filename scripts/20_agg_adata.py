from anndata import read_h5ad, AnnData
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
    parser.add_argument("--agg-func", type=str, required=True, choices=["mean", "sum"], default="sum", help = "Aggregation function to use (default: mean)")
    
    args = parser.parse_args()

    adata = read_h5ad(args.input_h5ad)

    adata.X = adata.layers["counts"] if adata.X is None else adata.X

    num_cells_orig = adata.shape[0]

    agg_func = args.agg_func

    unique_cell_types = adata.obs[args.cell_type_col].unique()
    unique_sample_ids = adata.obs[args.sample_id_col].unique()

    assert len(unique_cell_types) <= 1, f"Expected only one unique cell type in the input data, but found {len(unique_cell_types)}: {unique_cell_types}"
    assert len(unique_sample_ids) <= 1, f"Expected only one unique sample ID in the input data, but found {len(unique_sample_ids)}: {unique_sample_ids}"

    # Aggregate
    agg_X_1d = adata.X.sum(axis=0) if agg_func == "sum" else adata.X.mean(axis=0)
    # We need to expand dims to make it 2D again (1, num_genes) instead of 1D (num_genes,).
    agg_X = np.expand_dims(np.asarray(agg_X_1d).flatten(), axis=0)


    # For obs, we can aggregate any quantitative columns, and for string/categorical column: either take the first value (if all uniform) or specify a special "mixed values" string.
    def agg_column(x):
        if pd.api.types.is_numeric_dtype(x):
            return x.sum() if agg_func == "sum" else x.mean()
        else:
            return x.iloc[0] if x.nunique() == 1 else "mixed values"
    
    if adata.obs.shape[0] == 0:
        first_row = {
            col: np.nan for col in adata.obs.columns
        }
        first_row = {
            **first_row,
            args.cell_type_col: args.cell_type_name,
            args.sample_id_col: args.sample_id,
        }
        agg_obs = pd.DataFrame(columns=adata.obs.columns, data=[first_row])
    else:
        agg_obs = (
            adata.obs
                .groupby([args.cell_type_col, args.sample_id_col])
                .agg(agg_column)
                .reset_index()
        )

    # Append column for the number of cells that were aggregated, which can be used for filtering later.
    agg_obs["num_cells_orig"] = num_cells_orig

    print(agg_obs)

    # No need to do anything with var
    agg_var = adata.var.copy()

    # Create new AnnData object
    agg_adata = AnnData(X=agg_X, obs=agg_obs, var=agg_var)

    # Save
    agg_adata.write_h5ad(args.output_h5ad)


