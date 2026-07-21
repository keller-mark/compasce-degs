from compasce_degs.io import MultiComparisonMetadata, write_zdone
import zarr
import argparse
from filelock import SoftFileLock
import os
from os.path import join

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--zarr-path", type=str, required=True, help = "Path to zarr store.")
    args = parser.parse_args()

    zarr_path = args.zarr_path

    # Find all of the suffixes of the form {zarr_path}/uns/comparison_metadata.{suffix}
    uns_files = os.listdir(join(zarr_path, "uns"))
    metadata_suffixes = [
        filename[20:] for filename in uns_files
        if filename.startswith("comparison_metadata.") and filename != "comparison_metadata.merged"
    ]

    cm = MultiComparisonMetadata()
    cm.load_state(zarr_path, include_comparisons=True)
    cm.merge_states(zarr_path, suffixes=metadata_suffixes)

    zarr_dir_lock = SoftFileLock(f"{zarr_path}.lock", poll_interval=1.0)
    with zarr_dir_lock:
        z = zarr.open(zarr_path, mode="a")
        z["/uns/comparison_metadata"] = cm.serialize()
        z["/uns/comparison_metadata"].attrs["encoding-type"] = "string"
        z["/uns/comparison_metadata"].attrs["encoding-version"] = "0.2.0"
        write_zdone(zarr_path, arr_path=["uns", "comparison_metadata.merged"])