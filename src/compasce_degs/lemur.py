import scanpy as sc
import pandas as pd
import numpy as np
from anndata import AnnData
import pylemur


def compute_lemur(ladata, cm):

    if ladata.has_zdone(["uns", "compute_lemur"]):
        return ladata

    sample_group_pairs = cm.sample_group_pairs

    def get_input_arr():
        return ladata.get_da_from_zarr_layer("logcounts")
    
    # Dict of { col1: [(l1, r1), (l2, r2), ...] }
    sample_group_pairs_by_col = {}
    for sample_group_pair in sample_group_pairs:
        sample_group_col, (sample_group_left, sample_group_right) = sample_group_pair
        if sample_group_col in sample_group_pairs_by_col:
            sample_group_pairs_by_col[sample_group_col].append((sample_group_left, sample_group_right))
        else:
            sample_group_pairs_by_col[sample_group_col] = [(sample_group_left, sample_group_right)]


    for sample_group_col, sample_group_pair_list in sample_group_pairs_by_col.items():
        cmp = cm.add_comparison([("compare", sample_group_col)])
        model = pylemur.tl.LEMUR(ladata, get_input_arr, design = f"~ {sample_group_col}", n_embedding=15, layer = "logcounts", copy=False)
        # model = pylemur.tl.LEMUR(ladata, design = f"~ {sample_group_col}", n_embedding=15, layer = "logcounts", copy=False)
        model.fit()
        model.align_with_harmony()

        method_params = {
            "design": f"~ {sample_group_col}",
            "n_embedding": 15,
            "layer": "logcounts"
        }

        # Recalculate the DensMAP on the embedding calculated by LEMUR
        ladata.obsm["lemur_embedding"] = model.embedding
        ladata.save(arr_path=["obsm", "lemur_embedding"])

        # Map a particular layer to X, to trick the scanpy functions that only work with X.
        ladata.set_alias(on_disk_path=["obsm", "lemur_embedding"], in_mem_path=["X"])

        sc.pp.neighbors(ladata, n_neighbors=30, use_rep="X", key_added="lemur_embedding_neighbors")
        sc.tl.umap(ladata, method="densmap", key_added="lemur_densmap", neighbors_key="lemur_embedding_neighbors")

        ladata.save(arr_path=["obsm", "lemur_densmap"])

        obsm_key = cmp.append_df("obsm", "lemur_embedding", method_params, {
            "obsType": "cell",
            "featureType": "gene",
            "sampleSetGroup": [[sample_group_col]],
        })
        ladata.obsm[obsm_key] = model.embedding

        obsm_key = cmp.append_df("obsm", "lemur_densmap", method_params, {
            "obsType": "cell",
            "featureType": "gene",
            "sampleSetGroup": [[sample_group_col]],
        })
        ladata.obsm[obsm_key] = ladata.obsm["lemur_densmap"]
        
    
        for sample_group_pair in sample_group_pair_list:
            (sample_group_left, sample_group_right) = sample_group_pair
            cmp = cm.add_comparison([("compare", sample_group_col), ("val", sample_group_left), ("val", sample_group_right)])

            ctrl_pred = model.predict(new_condition=model.cond(**{ sample_group_col: sample_group_left }))
            stim_pred = model.predict(new_condition=model.cond(**{ sample_group_col: sample_group_right }))
            lemur_diff_matrix = (stim_pred - ctrl_pred)
            
            # Store the results
            layer_key = cmp.append_df("layers", "lemur_diff", method_params, {
                "obsType": "cell",
                "featureType": "gene",
                "sampleSetFilter": [[sample_group_col, sample_group_left], [sample_group_col, sample_group_right]],
            })
            ladata.layers[layer_key] = lemur_diff_matrix
    
    ladata.write_zdone(["uns", "compute_lemur"])

    return ladata
