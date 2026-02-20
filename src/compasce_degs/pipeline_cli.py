import argparse

from ._o2 import create_o2_dask_client
from .normalization import normalize_basic, normalize_pearson_residuals
from .densmap import densmap
from .diffexp import compute_diffexp
from .diffexp_pydeseq2 import compute_diffexp_pydeseq2
from .diffabundance import compute_diffabundance
from .lemur import compute_lemur
from .io.lazy_anndata import LazyAnnData
from .io.comparison_metadata import MultiComparisonMetadata


def run_cli():
    parser = argparse.ArgumentParser()
    parser.add_argument("--zarr-path", type=str, help="Path to zarr store")
    parser.add_argument("--input-deg-dir", type=str, required=False, help = "Path to folder containing precomputed DEG .txt files.")
    parser.add_argument("--mem-limit", type=str, default='16GB', required=False)
    parser.add_argument("--n-workers", type=int, default=2, required=False)
    parser.add_argument("--threads-per-worker", type=int, default=2, required=False)
    parser.add_argument("--function-name", type=str, choices=[
        'normalize_basic',
        'normalize_pearson_residuals',
        'densmap',
        'compute_diffexp',
        'compute_diffexp_pydeseq2',
        'compute_diffabundance',
        'compute_lemur',
    ])
    args = parser.parse_args()
    
    zarr_path = args.zarr_path
    function_name = args.function_name

    func_mapping = {
        'normalize_basic': normalize_basic,
        'normalize_pearson_residuals': normalize_pearson_residuals,
        'densmap': densmap,
        'compute_diffexp': compute_diffexp,
        'compute_diffexp_pydeseq2': compute_diffexp_pydeseq2,
        'compute_diffabundance': compute_diffabundance,
        'compute_lemur': compute_lemur,
    }

    func_to_run = func_mapping[function_name]
    client = create_o2_dask_client(
        memory_limit=args.mem_limit,
        n_workers=args.n_workers,
        threads_per_worker=args.threads_per_worker,
    )

    cm = MultiComparisonMetadata()
    cm.load_state(zarr_path, include_comparisons=False)

    ladata = LazyAnnData(zarr_path, client=client)

    kwargs = {}
    if function_name == "compute_diffexp":
        kwargs["input_deg_dir"] = args.input_deg_dir

    func_to_run(ladata, cm, **kwargs)

    ladata.uns[f"comparison_metadata.{function_name}"] = cm.serialize()
    ladata.save(arr_path=["uns", f"comparison_metadata.{function_name}"])

