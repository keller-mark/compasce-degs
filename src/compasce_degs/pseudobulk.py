import scanpy as sc
import numpy as np
import pandas as pd
import dask.array as da
from anndata import AnnData

from .dask import create_dask_wrapper


# Copied from https://github.com/scverse/decoupler/blob/b232862067217adbc87941a4337be80f1fa984f0/src/decoupler/pp/anndata.py#L164
def _sample_group_by(
    sample: str,
    group: str | None,
    obs: pd.DataFrame,
    verbose: bool = False,
):
    # Use one column if the same
    if sample == group:
        group = None
    # Handle list columns
    ocols = obs.select_dtypes(include="object").columns
    for ocol in ocols:
        has_list = any(isinstance(x, list) for x in obs[ocol].values)
        if has_list:
            obs[ocol] = obs[ocol].str.join("_")
    if group is None:
        # Filter extra columns in obs
        cols = obs.groupby(sample, observed=True).nunique(dropna=False).eq(other=1).all(axis=0)
        cols = np.hstack([sample, cols[cols].index])
        obs = obs.loc[:, cols]
        # Get unique samples
        smples = obs[sample].unique()
        groups = None
        # Get number of obs
        n_rows = len(smples)
        gsize = 0
    else:
        # Check if extra grouping is needed
        if type(group) is list:
            joined_cols = "_".join(group)
            obs[joined_cols] = obs[group[0]].str.cat(obs[group[1:]].astype("U"), sep="_")
            group = joined_cols
        # Filter extra columns in obs
        cols = obs.groupby([sample, group], observed=True).nunique(dropna=False).eq(other=1).all(axis=0)
        cols = np.hstack([sample, group, cols[cols].index])
        obs = obs.loc[:, cols]
        # Get unique samples and groups
        smples = np.unique(obs[sample].values)
        groups = np.unique(obs[group].values)
        # Get number of obs
        n_rows = len(smples) * len(groups)
        gsize = groups.size
    m = f"Generating {n_rows} profiles: {smples.size} samples x {gsize} groups"
    return obs, group, smples, groups, n_rows



def _pseudobulk(X, n_rows, n_cols, sample_col, groups_col, samples, groups, obs, new_obs, mode="sum"):
    # Reference: https://github.com/scverse/decoupler/blob/b232862067217adbc87941a4337be80f1fa984f0/src/decoupler/pp/anndata.py#L305

    # Init empty variables
    # psbulk = np.zeros((n_rows, n_cols))
    # props = np.zeros((n_rows, n_cols))
    ncells = np.zeros(n_rows)
    counts = np.zeros(n_rows)

    # Iterate for each group and sample
    i = 0
    psbulk_rows = []
    for grp in groups:
        for smp in samples:
            # Get cells from specific sample and group
            # Reference: https://github.com/scverse/decoupler/blob/b232862067217adbc87941a4337be80f1fa984f0/src/decoupler/pp/anndata.py#L264
            obs_mask = ((obs[sample_col] == smp) & (obs[groups_col] == grp)).values
            # Select rows and sum
            group_sum = da.sum(X[obs_mask, :], axis=0)

            # Skip if few cells or not enough counts
            ncell = group_sum.shape[0]
            count = da.sum(group_sum).compute()
            ncells[i] = ncell
            counts[i] = count
            # m = f"group={grp}\tsample={smp}\tcells={ncell}\tcounts={count}"
            # _log(m, level="info", verbose=verbose)
            # Write new meta-data
            index = smp + "_" + grp
            tmp = obs[(obs[sample_col] == smp) & (obs[groups_col] == grp)].drop_duplicates().values
            if tmp.shape[0] == 0:
                tmp = obs[obs[sample_col] == smp].drop(columns=groups_col).drop_duplicates()
                tmp = tmp.head(1)  # Remove extra repeated cat variables
                tmp[groups_col] = grp
                tmp = tmp[obs.columns].values
            new_obs.loc[index, :] = tmp
            psbulk_rows.append(group_sum)
            i += 1
    # Stack all rows together to form final dask array.
    final_dask_array = da.stack(psbulk_rows, axis=0)
    return final_dask_array


def pseudobulk(ladata, sample_col, groups_col, layer="counts", mode="sum"):
    var_df = ladata.var.copy() # unchanged in resulting anndata object
    obs_df = ladata.obs.copy()
    
    obs, groups_col, samples, groups, n_rows = _sample_group_by(
        sample=sample_col,
        group=groups_col,
        obs=obs_df,
    )
    n_cols = var_df.index.size
    new_obs = pd.DataFrame(columns=obs.columns)
    
    def get_input_arr():
        return ladata.get_da_from_zarr_layer(layer)
    
    def put_output_arr(output_arr):
        # Make cats
        for col in new_obs.columns:
            if not pd.api.types.is_numeric_dtype(new_obs[col]):
                new_obs[col] = new_obs[col].astype("category")

        X = output_arr.compute()

        print("Pseudobulked matrix shape:", X.shape)
        
        # Create new AnnData
        psbulk = AnnData(X=X, obs=new_obs, var=var_df)
        return psbulk

    pseudobulk_dask = create_dask_wrapper(_pseudobulk)
    return pseudobulk_dask(get_input_arr, put_output_arr, mode=mode, n_rows=n_rows, n_cols=n_cols, sample_col=sample_col, groups_col=groups_col, samples=samples, groups=groups, obs=obs, new_obs=new_obs)
