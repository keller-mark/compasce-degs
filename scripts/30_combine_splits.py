from anndata import read_h5ad, AnnData
import numpy as np
import pandas as pd

if __name__ == "__main__":
    # No argparse; snakemake script.

    # Get list of input files from snakemake
    input_files = snakemake.input
    output_file = snakemake.output[0]

    # Read all input files and concatenate their X and obs.
    # Just take the first var, assume they are all the same.

    adatas = [read_h5ad(f) for f in input_files]

    X_concat = np.concatenate([adata.X for adata in adatas], axis=0)
    obs_concat = pd.concat([adata.obs for adata in adatas], axis=0)
    var = adatas[0].var.copy()

    # Concatenate the layers as well.
    # Only keep layers that are present in all inputs.
    layer_keys = [k for k in adatas[0].layers.keys() if all(k in adata.layers for adata in adatas)]
    layers_concat = {
        k: np.concatenate([np.asarray(adata.layers[k]) for adata in adatas], axis=0)
        for k in layer_keys
    }

    # Create new AnnData object
    combined_adata = AnnData(X=X_concat, obs=obs_concat, var=var, layers=layers_concat)
    combined_adata.write_h5ad(output_file)
