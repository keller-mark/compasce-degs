import argparse
from os.path import basename, dirname

from ._o2 import create_o2_dask_client
from .normalization import normalize_basic, normalize_pearson_residuals
from .densmap import densmap
from .diffexp import compute_diffexp
from .diffexp_pydeseq2 import compute_diffexp_pydeseq2
from .diffabundance import compute_diffabundance
from .lemur import compute_lemur
from .insert_deg_results import insert_celltype_vs_rest_degs, insert_within_celltype_case_vs_control_degs
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

    # Add subparsers for commands that need their own arguments.
    subparsers = parser.add_subparsers(dest='command', help='subcommand help')

    # Create the parser for the "insert_celltype_vs_rest_degs" command
    subparser_insert_celltype_vs_rest_degs = subparsers.add_parser('insert_celltype_vs_rest_degs')
    subparser_insert_celltype_vs_rest_degs.add_argument("--csv-path", type=str, required=True, help = "Path to input CSV file.")
    subparser_insert_celltype_vs_rest_degs.add_argument("--cell-type-col", type=str, required=True, help = "Name of cell type column")
    subparser_insert_celltype_vs_rest_degs.add_argument("--sample-id-col", type=str, required=True, help = "Name of sample ID column")
    subparser_insert_celltype_vs_rest_degs.add_argument("--cell-type-name", type=str, required=True, help = "Cell type to subset for")
    subparser_insert_celltype_vs_rest_degs.add_argument("--agg-func", type=str, required=True, help = "Aggregation function used")
    subparser_insert_celltype_vs_rest_degs.add_argument("--out-path", type=str, required=True, help = "Output path")


    # Create the parser for the "insert_within_celltype_case_vs_control_degs" command
    subparser_within_celltype_degs = subparsers.add_parser('insert_within_celltype_case_vs_control_degs')
    subparser_within_celltype_degs.add_argument("--csv-path", type=str, required=True, help = "Path to input CSV file.")
    subparser_within_celltype_degs.add_argument("--cell-type-col", type=str, required=True, help = "Name of cell type column")
    subparser_within_celltype_degs.add_argument("--sample-id-col", type=str, required=True, help = "Name of sample ID column")
    subparser_within_celltype_degs.add_argument("--cell-type-name", type=str, required=True, help = "Cell type to subset for")
    subparser_within_celltype_degs.add_argument("--sample-group-col", type=str, required=True, help = "Name of sample group column")
    subparser_within_celltype_degs.add_argument("--sample-group-lhs", type=str, required=True, help = "Value of sample group LHS")
    subparser_within_celltype_degs.add_argument("--sample-group-rhs", type=str, required=True, help = "Value of sample group RHS")
    subparser_within_celltype_degs.add_argument("--agg-func", type=str, required=True, help = "Aggregation function used")
    subparser_within_celltype_degs.add_argument("--out-path", type=str, required=True, help = "Output path")

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
        'insert_celltype_vs_rest_degs': insert_celltype_vs_rest_degs,
        'insert_within_celltype_case_vs_control_degs': insert_within_celltype_case_vs_control_degs,
    }

    out_key = f"comparison_metadata.{function_name}"
    if args.function_name is None:
        out_key = basename(dirname(args.out_path)) # we want "something" from /path/to/something/.zdone
        
        if args.command == "insert_celltype_vs_rest_degs":
            function_name = "insert_celltype_vs_rest_degs"
        elif args.command == "insert_within_celltype_case_vs_control_degs":
            function_name = "insert_within_celltype_case_vs_control_degs"
        else:
            raise ValueError("Either --function-name or --command must be provided.")
    
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
    
    if args.command == "insert_celltype_vs_rest_degs":
        kwargs["csv_path"] = args.csv_path
        kwargs["cell_type_col"] = args.cell_type_col
        kwargs["sample_id_col"] = args.sample_id_col
        kwargs["cell_type_name"] = args.cell_type_name
        kwargs["agg_func"] = args.agg_func
        kwargs["out_key"] = out_key
    elif args.command == "insert_within_celltype_case_vs_control_degs":
        kwargs["csv_path"] = args.csv_path
        kwargs["cell_type_col"] = args.cell_type_col
        kwargs["sample_id_col"] = args.sample_id_col
        kwargs["cell_type_name"] = args.cell_type_name
        kwargs["sample_group_col"] = args.sample_group_col
        kwargs["sample_group_lhs"] = args.sample_group_lhs
        kwargs["sample_group_rhs"] = args.sample_group_rhs
        kwargs["agg_func"] = args.agg_func
        kwargs["out_key"] = out_key

    func_to_run(ladata, cm, **kwargs)

    ladata.uns[out_key] = cm.serialize()
    ladata.save(arr_path=["uns", out_key])

