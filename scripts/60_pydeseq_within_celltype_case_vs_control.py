from anndata import read_h5ad, AnnData
import numpy as np
import pandas as pd
import pertpy as pt
import argparse

NUM_CELLS_COLNAME = "num_cells_orig"

def clean_deg_df(df):
    # Rename the pds2.test_contrasts output to match scanpy's rank_genes_groups_df format
    # Scanpy colnames = ["names", "scores", "logfoldchanges", "pvals", "pvals_adj"]
    # References:
    # - https://pertpy.readthedocs.io/en/stable/tutorials/notebooks/differential_gene_expression.html#differential-expression-testing-with-pydeseq2
    # - https://scanpy.readthedocs.io/en/stable/generated/scanpy.get.rank_genes_groups_df.html
    df = df.rename(columns={
        "log_fc": "logfoldchanges",
        "p_value": "pvals",
        "adj_p_value": "pvals_adj",
        # TODO: rename the other columns as well?
        # baseMean
        # lfcSE
        # stat
    })
    df = df.sort_values(by="pvals_adj", ascending=True)
    df = df.set_index("variable")
    return df



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-h5ad", type=str, required=True, help = "Path to input H5AD file.")
    parser.add_argument("--output-de-csv", type=str, required=True, help = "Path to output CSV file for DE results.")
    parser.add_argument("--output-obs-filtering-csv", type=str, required=True, help = "Path to output CSV file for filtering obs results, for debugging.")
    parser.add_argument("--output-var-filtering-csv", type=str, required=True, help = "Path to output CSV file for filtering var results, for debugging.")
    parser.add_argument("--cell-type-col", type=str, required=True, help = "Name of cell type column")
    parser.add_argument("--cell-type-name", type=str, required=True, help = "Cell type to subset for")
    parser.add_argument("--sample-id-col", type=str, required=True, help = "Name of sample ID column")

    parser.add_argument("--sample-group-col", type=str, required=True, help = "Name of sample group column")
    parser.add_argument("--sample-group-lhs", type=str, required=True, help = "Left-hand side of sample group design (e.g. 'Primary Adjudicated Category')")
    parser.add_argument("--sample-group-rhs", type=str, required=True, help = "Right-hand side of sample group design (e.g. 'AKI')")

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

    sample_group_col = args.sample_group_col
    sample_group_lhs = args.sample_group_lhs
    sample_group_rhs = args.sample_group_rhs

    # Filter to the current cell type
    pdata = pdata[pdata.obs[cell_type_col] == cell_type].copy()

    # Filter to the pair of sample groups
    pdata = pdata[pdata.obs[sample_group_col].isin([sample_group_lhs, sample_group_rhs])].copy()


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

    pdata.obs["has_sufficient_num_samples"] = pdata.obs[cell_type_col].apply(lambda cell_type_val: num_samples_per_cell_type[cell_type_val] >= num_samples_threshold)

    pdata.obs["has_sufficient_num_samples"] = pdata.obs["has_sufficient_num_samples"].astype(bool)
    pdata.obs["has_sufficient_num_cells"] = pdata.obs["has_sufficient_num_cells"].astype(bool)

    # Save this filtering info and write for debugging
    # TODO: subset filtering df to only the current cell type? otherwise, all the output files are redundant.
    filtering_by_num_cells_df = pdata.obs[["has_sufficient_num_cells", "has_sufficient_num_samples", cell_type_col, sample_id_col, sample_group_col, NUM_CELLS_COLNAME]].copy()
    filtering_by_num_cells_df.to_csv(args.output_obs_filtering_csv, index=True)

    # Subset the anndata object.
    pdata = pdata[pdata.obs["has_sufficient_num_cells"] & pdata.obs["has_sufficient_num_samples"]].copy()

    # Filter along var: remove genes not expressed (count > 0) in at least 50% of pseudobulked samples.
    n_samples = pdata.X.shape[0]
    num_expressed = np.asarray((pdata.X > 0).sum(axis=0)).flatten()
    frac_expressed = num_expressed / n_samples
    var_mask = frac_expressed >= 0.5

    # TODO: also manually (via regex) remove "AC" and "AL"-prefixed genes?

    filtering_by_var_df = pdata.var.copy()
    filtering_by_var_df["frac_expressed"] = frac_expressed
    filtering_by_var_df["is_expressed_in_sufficient_samples"] = var_mask
    filtering_by_var_df.to_csv(args.output_var_filtering_csv, index=True)

    pdata = pdata[:, var_mask].copy()

    print(f"Filtered pdata shape: {pdata.shape}")

    # LHS vs RHS
    has_lhs_and_rhs = pdata.obs[sample_group_col].nunique() == 2
    if not has_lhs_and_rhs:
        print(f"Warning: after filtering, there are not at least 2 groups to compare for cell type {cell_type} and {sample_group_col}: {sample_group_lhs} vs {sample_group_rhs}. Skipping DE analysis for this pair and cell type.")
        empty_df = pd.DataFrame(columns=["variable", "logfoldchanges", "pvals", "pvals_adj"])
        empty_df.to_csv(args.output_de_csv, index=False)
    else:
        # For cell type vs. rest, we use a design such as "~subclass_l1"
        # Reference: https://hbctraining.github.io/DGE_workshop/lessons/04_DGE_DESeq2_analysis.html
        pds2 = pt.tl.PyDESeq2(adata=pdata, design=f"~{sample_group_col}")
        pds2.fit()

        df = pds2.test_contrasts(pds2.contrast(column=sample_group_col, baseline=sample_group_rhs, group_to_compare=sample_group_lhs))
        df = clean_deg_df(df)

        print("Done with DE analysis, writing output CSV...")

        df.to_csv(args.output_de_csv, index=True)
