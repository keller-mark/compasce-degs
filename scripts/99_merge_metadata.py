from compasce_degs.io import MultiComparisonMetadata, write_zdone
import zarr
import argparse
from filelock import SoftFileLock
import os
from os.path import join

# Result types written to /uns/comparison_metadata by the initial (clean_h5ad) step.
# Everything else comes from the per-step /uns/comparison_metadata.{suffix} entries.
BASE_COMPARISON_KEY = "__all__"
BASE_ANALYSIS_TYPES = {"cells", "samples"}


def get_base_comparisons(prev_comparisons):
    # /uns/comparison_metadata may already contain the output of a previous merge.
    # Since deep_update concatenates the "results" lists, merging on top of a previous
    # merge would duplicate entries and keep stale entries whose underlying data has been
    # deleted (e.g., to be re-computed). To keep this script idempotent, start from only
    # the base results that were written by the initial step.
    if BASE_COMPARISON_KEY not in prev_comparisons:
        raise ValueError(f"Expected '{BASE_COMPARISON_KEY}' in /uns/comparison_metadata comparisons.")
    base = prev_comparisons[BASE_COMPARISON_KEY]
    base_results = []
    seen_paths = set()
    for result in base["results"]:
        if result["analysis_type"] in BASE_ANALYSIS_TYPES and result["path"] not in seen_paths:
            base_results.append(result)
            seen_paths.add(result["path"])
    return {
        BASE_COMPARISON_KEY: {
            **base,
            "results": base_results,
        }
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--zarr-path", type=str, required=True, help = "Path to zarr store.")
    args = parser.parse_args()

    zarr_path = args.zarr_path

    # Find all of the suffixes of the form {zarr_path}/uns/comparison_metadata.{suffix}
    uns_files = os.listdir(join(zarr_path, "uns"))
    metadata_suffixes = sorted([
        filename[len("comparison_metadata."):] for filename in uns_files
        if filename.startswith("comparison_metadata.") and filename != "comparison_metadata.merged"
    ])

    cm = MultiComparisonMetadata()
    cm.load_state(zarr_path, include_comparisons=True)
    cm._prev_comparisons_dict = get_base_comparisons(cm._prev_comparisons_dict)
    cm.merge_states(zarr_path, suffixes=metadata_suffixes)

    zarr_dir_lock = SoftFileLock(f"{zarr_path}.lock", poll_interval=1.0)
    with zarr_dir_lock:
        z = zarr.open(zarr_path, mode="a")
        # Note: this replaces the array directory, which removes the .zdone file
        # inside it (an output of the clean_h5ad Snakemake rule), so we re-write it below.
        z["/uns/comparison_metadata"] = cm.serialize()
        z["/uns/comparison_metadata"].attrs["encoding-type"] = "string"
        z["/uns/comparison_metadata"].attrs["encoding-version"] = "0.2.0"
        write_zdone(zarr_path, arr_path=["uns", "comparison_metadata"])
        write_zdone(zarr_path, arr_path=["uns", "comparison_metadata.merged"])
