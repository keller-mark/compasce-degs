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
    print(f"Inserting pseudobulk for cell types vs rest, using PyDESeq2")

    # Check for a .zdone file
    if ladata.has_zdone(["uns", out_key]):
        return ladata
    
    # Check whether df is empty (insufficient data for DEG analysis).
    df = pd.read_csv(csv_path)
    is_empty = df.shape[0] == 0

    cmp = cm.add_comparison([("compare", cell_type_col), ("val", cell_type_name), "__rest__"])

    uns_key = cmp.append_df("uns", "pydeseq2", {
        "contrast": {
            "column": "is_cell_type",
            "baseline": "rest",
            "group_to_compare": cell_type_name,
        },
        "is_empty": is_empty,
    }, {
        "obsType": "cell",
        "featureType": "gene",
        "obsSetSelection": [[cell_type_col, cell_type_name]],
    })

    if not is_empty:
        ladata.uns[uns_key] = df

    # Write a .zdone file
    ladata.write_zdone(["uns", out_key])

    return ladata






def insert_within_celltype_case_vs_control_degs(
    ladata, cm,
    csv_path, cell_type_col, sample_id_col, cell_type_name, sample_group_col, sample_group_lhs, sample_group_rhs, agg_func,
    out_key,
):
    print(f"Inserting pseudobulk for within cell type, case vs control, using PyDESeq2")

    # Check for a .zdone file
    if ladata.has_zdone(["uns", out_key]):
        return ladata
    
    # Check whether df is empty (insufficient data for DEG analysis).
    df = pd.read_csv(csv_path)
    is_empty = df.shape[0] == 0

    cmp = cm.add_comparison([("filter", cell_type_col), ("val", cell_type_name), ("compare", sample_group_col), ("val", sample_group_lhs), ("val", sample_group_rhs)])

    uns_key = cmp.append_df("uns", "pydeseq2", {
        "contrast": {
            "column": "cell_type_sample_group",
            "baseline": f"{cell_type_name}_{sample_group_lhs}",
            "group_to_compare": f"{cell_type_name}_{sample_group_rhs}",
        },
        "is_empty": is_empty,
    }, {
        "obsType": "cell",
        "featureType": "gene",
        "obsSetFilter": [[cell_type_col, cell_type_name]],
        "sampleSetSelection": [[sample_group_col, sample_group_rhs]],
        "sampleSetFilter": [[sample_group_col, sample_group_lhs], [sample_group_col, sample_group_rhs]],
    })

    if not is_empty:
        ladata.uns[uns_key] = df

    # Write a .zdone file
    ladata.write_zdone(["uns", out_key])

    return ladata