import os
from os.path import join
import zarr
import dask.array as da
from anndata import AnnData
from dask.distributed import progress

from .zarr_io import dispatched_read_zarr, dispatched_write_zarr, write_zdone as io_write_zdone, has_zdone as io_has_zdone, try_cast_arr


class MappingWrapper:
    # Because each value of `mapping` is a ZarrArray, we need to ensure that
    # writes are written to disk (instead of just updating the dictionary value for the key).

    def __init__(self, mapping, adata, parent_key):
        self.mapping = mapping
        self.adata = adata
        self.parent_key = parent_key

    def __getitem__(self, key):
        # TODO: return zarr array that 
        return self.mapping[key]

    def __contains__(self, key):
        return key in self.mapping
    
    def __setitem__(self, key, value):
        # assumes 2D
        if key in self.mapping:
            assert isinstance(self.mapping[key], zarr.Array)
            self.mapping[key][:, :] = value
        else:
            print(f"Setting {self.parent_key}/{key} via MappingWrapper. This will not be written to disk until ladata.save().")
            getattr(self.adata, self.parent_key)[key] = value
            # Do not write. Instead, let the user call ladata.save() to write to disk.
    
    def __getattr__(self, key):
        return getattr(self.mapping, key)

    def __delitem__(self, key):
        del self.mapping[key]

class LazyAnnData(AnnData):

    zarr_path = None
    var_chunk_size = None
    client = None
    
    adata = None
    z = None
    aliases = {}
    
    done_init = False

    def __init__(self, zarr_path, var_chunk_size=5, client=None):
        self.zarr_path = zarr_path
        self.var_chunk_size = var_chunk_size
        self.client = client

        self.adata = dispatched_read_zarr(zarr_path) # only contains var, obs, uns
        self.z = zarr.open(zarr_path, mode="a")
       
        super().__init__()

        self.done_init = True
    
    def has_zdone(self, arr_path):
        return io_has_zdone(self.zarr_path, arr_path=arr_path)

    def write_zdone(self, arr_path):
        io_write_zdone(self.zarr_path, arr_path=arr_path)


    def set_alias(self, on_disk_path, in_mem_path):
        # Set up an alias so that we can, for example, trick scanpy into
        # using a layer like z["/layers/counts"] for adata.X
        self.aliases[tuple(in_mem_path)] = list(on_disk_path)
    
    def clear_aliases(self):
        self.aliases = {}

    def __getattribute__(self, key):
        orig_getattr = object.__getattribute__
        done_init = orig_getattr(self, "done_init")
        if not done_init:
            return orig_getattr(self, key)
        
        z = orig_getattr(self, "z")
        # __getattr__ only gets called for attributes that don't actually exist
        #print(f"Getting {key}")
        if key in { "X" }:
            if ("X",) in self.aliases:
                on_disk_path = self.aliases[("X",)]
                return z[f"/{'/'.join(on_disk_path)}"]
            return z["/X"]
        elif key in { "obsm", "varm", "layers" }:
            on_disk_subkeys = z[key].keys()
            in_mem_subkeys = [tup[1] for tup in self.aliases if tup[0] == key]
            subkeys = set(on_disk_subkeys).union(in_mem_subkeys)

            print(f"Getting {key} subkeys: {subkeys}")
            def get_on_disk_path(subkey):
                if subkey in in_mem_subkeys:
                    on_disk_path = self.aliases[(key, subkey)]
                    return f"/{'/'.join(on_disk_path)}"
                return f"/{key}/{subkey}"

            return MappingWrapper({
                subkey: z[get_on_disk_path(subkey)]
                for subkey in subkeys
            }, self.adata, key)
        elif key in {"obs", "var", "uns", "obsp"}:
            return getattr(self.adata, key)
        
        return orig_getattr(self, key)
    
    def __setattr__(self, key, value):
        if not self.done_init:
            return super().__setattr__(key, value)
        print(f"Setting {key}")
        if key in {"obs", "var", "uns", "obsp"}:
            setattr(self.adata, key, value)
        elif key in { "obsm", "varm", "layers", "X" }:
            raise ValueError(f"Cannot set {key} via LazyAnnData. Set on the zarr array directly.")
        
        super().__setattr__(key, value)

    def save(self, arr_path=None, mode="a"):
        dispatched_write_zarr(self.adata, self.zarr_path, var_chunk_size=self.var_chunk_size, arr_path=arr_path, mode=mode, client=self.client)

    def copy_layer(self, src_key, dest_key):
        zarr.copy(self.z[f"/layers/{src_key}"], self.z["/layers"], name=dest_key)

    def get_da_from_zarr_layer(self, layer_key):
        X_path = f"/layers/{layer_key}"

        if self.client is None:
            raise ValueError("A Dask client must be provided to use get_da_from_zarr_layer")
            # TODO: instead of throwing, just run da.from_array() on the numpy array normally?

        z = zarr.open(self.zarr_path, path=X_path, mode='r')

        # Use Zarr so that X does not get loaded into memory
        # which would cause Dask to include it in the task graph,
        # which causes the task graph to be too large.
        X_dask = da.from_zarr(url=self.zarr_path, component=X_path, chunks=z.chunks, inline_array=False)
        return X_dask

    def put_da_to_zarr_layer(self, layer_key, X_dask):
        X_path = f"/layers/{layer_key}"

        if self.client is None:
            raise ValueError("A Dask client must be provided to use put_da_to_zarr_layer")
            # TODO: instead of throwing, just run .compute() and write the numpy array normally?

        # Store using Zarr in a distributed way so that we don't need to transfer data back from each worker,
        # which is too large for a single worker to handle.
        residuals_delayed = try_cast_arr(X_dask).to_zarr(url=self.zarr_path, component=X_path, overwrite=True, compute=False)
        residuals_future = self.client.persist(residuals_delayed)

        progress(residuals_future, notebook=False)
        #future_result = residuals_future.compute()

        self.write_zdone(["layers", layer_key])


def create_lazy_anndata(adata, zarr_path, client=None, overwrite=False):
    assert "counts" in adata.layers, "The AnnData object must have a 'counts' layer."

    dispatched_write_zarr(adata, zarr_path, arr_path=["layers", "counts"], mode=("w" if overwrite else "a"), client=client)

    return LazyAnnData(zarr_path, client=client)


def create_sample_df(adata, cm):
    sample_id_col = cm.sample_id_col
    assert sample_id_col in adata.obs.columns
    sample_groupby = adata.obs.groupby(by=sample_id_col)
    per_sample_obs_cols_nunique = sample_groupby.nunique()
    # Find columns whose values are only at most one unique value per sample.
    sample_cols = []
    for col in per_sample_obs_cols_nunique.columns:
        if len(per_sample_obs_cols_nunique[col].unique()) == 1 and per_sample_obs_cols_nunique[col].unique()[0] == 1:
            sample_cols.append(col)
    
    sample_df = sample_groupby.first()[sample_cols]
    return sample_df