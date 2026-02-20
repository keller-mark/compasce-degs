from os.path import join
import scanpy as sc
import pandas as pd
import numpy as np
from anndata import AnnData
from requests.exceptions import ConnectionError

import decoupler as dc
import pertpy as pt

from .pseudobulk import pseudobulk



# Functions for cleaning up dataframes
def cleanup_rank_genes_groups_df(df):
    # Rename the pds2.test_contrasts output to match scanpy's rank_genes_groups_df format
    # Scanpy colnames = ["names", "scores", "logfoldchanges", "pvals", "pvals_adj"]
    # References:
    # - https://pertpy.readthedocs.io/en/stable/tutorials/notebooks/differential_gene_expression.html#differential-expression-testing-with-pydeseq2
    # - https://scanpy.readthedocs.io/en/stable/generated/scanpy.get.rank_genes_groups_df.html
    df = df.rename(columns={
        "log2FoldChange": "logfoldchanges",
        "pvalue": "pvals",
        "padj": "pvals_adj",
    })
    df = df.sort_values(by="pvals_adj", ascending=True)
    df = df.set_index("names")
    return df


def compute_diffexp_pydeseq2(ladata, cm):
    print(f"Running pseudobulk for cell types vs rest, using PyDESeq2")

    # Check for a .zdone file
    if ladata.has_zdone(["uns", "compute_diffexp_pydeseq2"]):
        return ladata

    cell_type_cols = cm.cell_type_cols

    for cell_type_col in cell_type_cols:

        ladata.obs[cell_type_col] = ladata.obs[cell_type_col].astype(str).astype("category")

        cell_types = ladata.obs[cell_type_col].unique().tolist()
        cell_types = [x for x in cell_types if pd.notna(x)]
        
        # Reference: https://pertpy.readthedocs.io/en/stable/tutorials/notebooks/differential_gene_expression.html#pseudobulks
        # pdata = dc.pp.pseudobulk(ladata, sample_col="specimen", groups_col="subclass_l1", layer="counts", mode="sum", verbose=True)
        # pdata = dc.pp.pseudobulk(ladata, sample_col=cm.sample_id_col, groups_col=cell_type_col, layer="counts", mode="sum", empty=True, verbose=True)
        pdata = pseudobulk(ladata, sample_col=cm.sample_id_col, groups_col=cell_type_col, layer="counts", mode="sum")

        for cell_type in cell_types:
            print(f"Getting diffexp test results for {cell_type} vs rest")
            cmp = cm.add_comparison([("compare", cell_type_col), ("val", cell_type), "__rest__"])
            
            # In order to compare one cell type vs. the rest, we need to create a new column in pdata.obs to store this information
            pdata.obs["is_cell_type"] = pdata.obs[cell_type_col].apply(lambda x: cell_type if x == cell_type else "rest")

            # For cell type vs. rest, we use a design such as "~subclass_l1"
            # Reference: https://hbctraining.github.io/DGE_workshop/lessons/04_DGE_DESeq2_analysis.html
            pds2 = pt.tl.PyDESeq2(adata=pdata, design=f"~is_cell_type")
            pds2.fit()

            df = pds2.test_contrasts(pds2.contrast(column="is_cell_type", baseline="rest", group_to_compare=cell_type))

            #df = sc.get.rank_genes_groups_df(ladata, group=cell_type, key=key_added)
            df = cleanup_rank_genes_groups_df(df)

            uns_key = cmp.append_df("uns", "pydeseq2", {
                "contrast": {
                    "column": "is_cell_type",
                    "baseline": "rest",
                    "group_to_compare": cell_type,
                },
            }, {
                "obsType": "cell",
                "featureType": "gene",
                "obsSetSelection": [[cell_type_col, cell_type]],
            })
            ladata.uns[uns_key] = df

            print(df.head())

            # Clean up to reduce memory usage
            del pds2
            del df
        del pdata

        # Within cell type (case vs. control)
        sample_group_pairs = cm.sample_group_pairs
        
        for cell_type in cell_types:
            print(f"Running diffexp test for {cell_type} and sample group pairs")
            for sample_group_pair in sample_group_pairs:
                
                sample_group_col, (sample_group_left, sample_group_right) = sample_group_pair
                cmp = cm.add_comparison([("filter", cell_type_col), ("val", cell_type), ("compare", sample_group_col), ("val", sample_group_left), ("val", sample_group_right)])
                try:
                    ladata.obs["cell_type_sample_group"] = ladata.obs[cell_type_col].astype(str) + "_" + ladata.obs[sample_group_col].astype(str)
                    ladata.obs["cell_type_sample_group"] = ladata.obs["cell_type_sample_group"].astype(str).astype("category")
                    #sc.tl.rank_genes_groups(ladata, groupby="cell_type_sample_group", groups=[f"{cell_type}_{sample_group_right}"], reference=f"{cell_type}_{sample_group_left}", method="wilcoxon", layer="logcounts", key_added=key_added)
                    
                    # TODO: only pseudobulk once per unique cell_type_col and sample_group_col combination?
                    # pdata = dc.pp.pseudobulk(ladata, sample_col=cm.sample_id_col, groups_col="cell_type_sample_group", layer="counts", mode="sum", empty=True, verbose=True)
                    pdata = pseudobulk(ladata, sample_col=cm.sample_id_col, groups_col="cell_type_sample_group", layer="counts", mode="sum")
                    
                    # For cell type vs. rest, we use a design such as "~subclass_l1"
                    # Reference: https://hbctraining.github.io/DGE_workshop/lessons/04_DGE_DESeq2_analysis.html
                    pds2 = pt.tl.PyDESeq2(adata=pdata, design=f"~cell_type_sample_group")
                    pds2.fit()
                    
                    df = pds2.test_contrasts(pds2.contrast(column="cell_type_sample_group", baseline="{cell_type}_{sample_group_left}", group_to_compare="{cell_type}_{sample_group_right}"))
                    # df = sc.get.rank_genes_groups_df(ladata, group=f"{cell_type}_{sample_group_right}", key=key_added)
                    df = cleanup_rank_genes_groups_df(df)

                    uns_key = cmp.append_df("uns", "pydeseq2", {
                        "contrast": {
                            "column": "cell_type_sample_group",
                            "baseline": f"{cell_type}_{sample_group_left}",
                            "group_to_compare": f"{cell_type}_{sample_group_right}",
                        },
                    }, {
                        "obsType": "cell",
                        "featureType": "gene",
                        "obsSetFilter": [[cell_type_col, cell_type]],
                        "sampleSetSelection": [[sample_group_col, sample_group_right]],
                        "sampleSetFilter": [[sample_group_col, sample_group_left], [sample_group_col, sample_group_right]],
                    })
                    ladata.uns[uns_key] = df

                    del pdata
                    del pds2
                    del df

                except (IndexError, ValueError) as e:
                    print(f"Error: likely due to insufficient data for comparison for {cell_type} and sample group pair {sample_group_pair}")

        # TODO: within cell type (inside spatial region vs. outside)

    # Write a .zdone file
    ladata.write_zdone(["uns", "compute_diffexp_pydeseq2"])

    return ladata
