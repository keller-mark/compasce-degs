import numpy as np
import zarr
import json
import os
from os.path import join
import dask.array as da
from anndata import AnnData
from anndata._io.zarr import read_dataframe, _read_legacy_raw, _clean_uns
from anndata.experimental import read_dispatched, write_dispatched, read_elem
from dask.distributed import progress

def try_cast_arr(arr):
    if isinstance(arr, np.ndarray) or isinstance(arr, da.Array) or isinstance(arr, zarr.Array):
        if arr.dtype.kind in ['f', 'u', 'i'] and arr.dtype.itemsize == 8:
            # Cast 64-bit to 32-bit
            return arr.astype(f'<{arr.dtype.kind}4')
    return arr

def dispatched_read_zarr(store):
    # Function that reads an AnnData object from a Zarr store but omits certain keys.
    # Adapted from https://github.com/scverse/anndata/blob/1461fecd1712eefb1e5a5c0a75547b0e169a23d5/src/anndata/_io/zarr.py#L51
    if isinstance(store, zarr.Group):
        f = store
    else:
        f = zarr.open(store, mode="r")

    # Read with handling for backwards compat
    def callback(func, elem_name: str, elem, iospec):
        #print(f"Reading {elem_name}")
        if elem_name == "/X":
            return None
        if elem_name == "/raw":
            return None

        if elem_name == "/layers" or elem_name == "/obsm":
            # We want to trick the anndata read_basic_zarr function into thinking that this Zarr group
            # only contains a subset of keys.
            # Reference: https://github.com/scverse/anndata/blob/1461fecd1712eefb1e5a5c0a75547b0e169a23d5/src/anndata/_io/specs/methods.py#L145
            elem = {}

        if iospec.encoding_type == "anndata" or elem_name.endswith("/"):
            return AnnData(
                **{
                    k: read_dispatched(v, callback)
                    for k, v in elem.items()
                    if not k.startswith("raw.")
                }
            )
        elif elem_name.startswith("/raw."):
            return None
        elif elem_name in {"/obs", "/var"}:
            return read_dataframe(elem)
        #elif elem_name == "/raw":
            # Backwards compat
        #    return _read_legacy_raw(f, func(elem), read_dataframe, func)
        
        return func(elem)

    adata = read_dispatched(f, callback=callback)

    # Backwards compat (should figure out which version)
    if "raw.X" in f:
        raw = AnnData(**_read_legacy_raw(f, adata.raw, read_dataframe, read_elem))
        raw.obs_names = adata.obs_names
        adata.raw = raw

    # Backwards compat for <0.7
    if isinstance(f["obs"], zarr.Array):
        _clean_uns(adata)

    return adata

def has_zdone(out_path, arr_path=None):
    if os.path.isfile(join(out_path, *arr_path, ".zdone")):
        return True
    return False

def write_zdone(out_path, arr_path=None):
    # Write a hidden file at `out_path/arr_path/.zdone``
    # to indicate done-ness of an operation for Snakemake
    if arr_path is not None:
        os.makedirs(join(out_path, *arr_path), exist_ok=True)
        json_path = join(out_path, *arr_path, ".zdone")
        with open(json_path, 'w') as f:
            json.dump({"done": True}, f)

def dispatched_write_zarr(adata, out_path, var_chunk_size=5, arr_path=None, mode="a", client=None):
    # Write to Zarr and set custom chunk shape for layers
    # Reference: https://anndata.readthedocs.io/en/latest/tutorials/notebooks/%7Bread%2Cwrite%7D_dispatched.html
    def write_chunked(func, store, k, elem, dataset_kwargs, iospec):
        """Write callback that chunks X and layers"""

        def set_chunks(d, chunks=None):
            """Helper function for setting dataset_kwargs. Makes a copy of d."""
            d = dict(d)
            if chunks is not None:
                d["chunks"] = chunks
            else:
                d.pop("chunks", None)
            return d

        if iospec.encoding_type == "array":
            if 'layers' in k or k.endswith('X'):
                dataset_kwargs = set_chunks(dataset_kwargs, (adata.shape[0], var_chunk_size))
            else:
                dataset_kwargs = set_chunks(dataset_kwargs, None)

        if isinstance(elem, da.Array):
            raise ValueError("Dask arrays should not be passed to write_chunked")
        elif elem is None:
            print("Skipping writing of None element")
        else:
            elem = try_cast_arr(elem)
            try:
                func(store, k, elem, dataset_kwargs=dataset_kwargs)
            except zarr.errors.ContainsArrayError:
                print(f"zarr.errors.ContainsArrayError at {k}")

    z = zarr.open(out_path, mode=mode)

    # Monkey patch the clear method to prevent clearing the root group
    # Reference: https://github.com/scverse/anndata/blob/1461fecd1712eefb1e5a5c0a75547b0e169a23d5/src/anndata/_io/specs/registry.py#L299C1-L299C26
    old_clear = z.clear
    z.clear = (lambda: None) # Do not allow clearing the root group

    old_delitem = z.__class__.__delitem__
    def patched_delitem(self, item):
        if item == "/layers" or item == "/obsm" or item == "/uns" or item.startswith("/uns/comparison_metadata"):
            pass
        else:
            print(f"Deleting {item}")
            old_delitem(self, item)
    z.__class__.__delitem__ = patched_delitem

    write_dispatched(z, "/", adata, callback=write_chunked)
    # Restore (though not really necessary)
    z.clear = old_clear
    z.__class__.__delitem__ = old_delitem

    write_zdone(out_path, arr_path=arr_path)

