from .normalization import normalize_basic, normalize_pearson_residuals
from .densmap import densmap
from .diffexp import compute_diffexp
from .diffexp_pydeseq2 import compute_diffexp_pydeseq2
from .diffabundance import compute_diffabundance
from .lemur import compute_lemur
from .io.lazy_anndata import create_lazy_anndata, create_sample_df
from .io.comparison_metadata import MultiComparisonMetadata


def run_all(get_adata, out_h5ad_path, out_zarr_path, overwrite=False, client=None, sample_id_col=None, donor_id_col=None, sample_group_pairs=None, cell_type_cols=None, stop_early=False, input_deg_dir=None):
    """
    def get_adata():
        return read_h5ad("path/to/adata.h5ad")
    zarr_path = "path/to/adata.zarr"
    client = create_dask_client()

    run_all(get_adata, zarr_path, client=client)
    """

    # We ask the caller to provide a function which returns `adata`, 
    # so that we can free the memory within the function.
    adata = get_adata()
    
    cm = MultiComparisonMetadata(
        sample_group_pairs=sample_group_pairs,
        sample_id_col=sample_id_col,
        donor_id_col=donor_id_col,
        cell_type_cols=cell_type_cols if cell_type_cols is not None else ["cell_type"],
    )
    # if not overwrite:
    #     # Initialize with contents of /uns/comparison_metadata
    #     cm.load_state(zarr_path)

    all_cmp = cm.add_comparison("__all__")
    all_cmp.append_df("uns", "cells", None, { "obsType": "cell" }) # assumed

    if "counts" not in adata.layers:
        raise ValueError("adata.layers['counts'] must exist")
    if adata.shape[0] < 1:
        raise ValueError("adata must have at least one row")
    if adata.shape[1] < 1:
        raise ValueError("adata must have at least one column")

    # ladata = create_lazy_anndata(adata, zarr_path, client=client, overwrite=overwrite)

    #del adata

    if stop_early:
        adata.uns["comparison_metadata"] = cm.serialize()
        # ladata.save(arr_path=["uns", "comparison_metadata"])

        adata.write_h5ad(out_h5ad_path)

        # Write to zarr
        # Convert counts to dense
        adata.layers["counts"] = adata.layers["counts"].toarray()
        ladata = create_lazy_anndata(adata, out_zarr_path, client=client, overwrite=overwrite)

        del adata

        sample_df = create_sample_df(ladata, cm)
        uns_key = all_cmp.append_df("uns", "samples", None, {
            "obsType": "sample"
        })
        ladata.uns[uns_key] = sample_df

        ladata.uns["comparison_metadata"] = cm.serialize()
        ladata.save(arr_path=["uns", "comparison_metadata"])

        return ladata
    
    # TODO: for KPMP, fill in pre-processed differential expression results here? Or during the compute_diffexp step?
    # Or, just fill in at the end, after the normal pipeline has finished?

    # depends on: uns/write_metadata/layers/counts
    # creates: uns/write_metadata/layers/logcounts
    normalize_basic(adata, cm)

    adata.write_h5ad(out_h5ad_path)
    """ 
    #del ladata.layers["counts"]

    # depends on: uns/write_metadata/layers/counts
    # creates: /layers/pearson_residuals
    #normalize_pearson_residuals(ladata)

    del ladata.layers["logcounts"]

    densmap(ladata, cm)

    compute_diffexp(ladata, cm, input_deg_dir=input_deg_dir)

    ladata.uns["comparison_metadata"] = cm.serialize()
    ladata.save()

    compute_diffexp_pydeseq2(ladata, cm, input_deg_dir=input_deg_dir)

    ladata.uns["comparison_metadata"] = cm.serialize()
    ladata.save()

    compute_diffabundance(ladata, cm)

    ladata.uns["comparison_metadata"] = cm.serialize()
    ladata.save()

    print("Starting normalize_pearson_residuals")
    normalize_pearson_residuals(ladata, cm)
    #ladata.save()

    #print("Starting compute_lemur")
    #compute_lemur(ladata, cm)

    ladata.uns["comparison_metadata"] = cm.serialize()
    ladata.save()
    """

    return True