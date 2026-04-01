from os.path import join
import scanpy as sc
import pandas as pd
import numpy as np
from anndata import AnnData


def insert_celltype_vs_rest_degs(
    ladata, cm,
    csv_path, cell_type_col, sample_id_col, cell_type_name, agg_func,
    out_key,
):
    print(f"Running pseudobulk for cell types vs rest, using PyDESeq2")

    # Check for a .zdone file
    if ladata.has_zdone(["uns", out_key]):
        return ladata

    cmp = cm.add_comparison([("compare", cell_type_col), ("val", cell_type_name), "__rest__"])

    uns_key = cmp.append_df("uns", "pydeseq2", {
        "contrast": {
            "column": "is_cell_type",
            "baseline": "rest",
            "group_to_compare": cell_type_name,
        },
    }, {
        "obsType": "cell",
        "featureType": "gene",
        "obsSetSelection": [[cell_type_col, cell_type_name]],
    })

    df = pd.read_csv(csv_path)

    ladata.uns[uns_key] = df