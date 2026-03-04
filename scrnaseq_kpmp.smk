include: "./common.smk"
configfile: "./scrnaseq_kpmp.yaml"

RAW_H5AD_PATH = join(RAW_DIR, "kpmp-aug-2025", "SingleNucleus_KPMP_Explorer_05182025.h5ad")
RAW_SAMPLES_PATH = join(RAW_DIR, "kpmp-aug-2025", "20250606_OpenAccessClinicalData.csv")

# Intermediate output paths
CLEANED_H5AD_PATH = join(INTERMEDIATE_DIR, "cleaned.h5ad")
NORMALIZED_H5AD_PATH = join(INTERMEDIATE_DIR, "normalized.h5ad")

L1_CELL_TYPES = sorted(config["cell_types"]["subclass_l1"])
L2_CELL_TYPES = sorted(config["cell_types"]["subclass_l2"])
L3_CELL_TYPES = sorted(config["cell_types"]["subclass_l3"])

def cell_type_name_to_index(cell_type_name, cell_type_col):
  # Convert cell type name to index in config["cell_types"]["subclass_l1"], which
  # is the order we use for the cell type columns in the .adata.h5ad files.
  return config["cell_types"][cell_type_col].index(cell_type_name)

def cell_type_index_to_name(cell_type_index, cell_type_col):
  return config["cell_types"][cell_type_col][cell_type_index]
  
def normalize_cell_type(cell_type):
    # Normalize cell type names to be filesystem-friendly, since we'll be using them in file paths.
    return (
      cell_type
        .replace("/", "__slash__")
        .replace(" ", "__space__")
        .replace("+", "__plus__")
        .replace("\xEF", "__i__")
    )

def unnormalize_cell_type(cell_type):
    return (
      cell_type
        .replace("__slash__", "/")
        .replace("__space__", " ")
        .replace("__plus__", "+")
        .replace("__i__", "\xEF")
    )

SPECIMEN_ID_COL = "specimen"
SPECIMEN_IDS = sorted(config["specimen_ids"])

def specimen_id_to_index(specimen_id):
  return SPECIMEN_IDS.index(specimen_id)

def index_to_specimen_id(specimen_id_index):
  return SPECIMEN_IDS[specimen_id_index]

NUM_SAMPLES_THRESHOLD = config["thresholds"]["min_samples_per_group"] # Min number of samples when using pseudobulked data
NUM_CELLS_PER_SAMPLE_THRESHOLD = config["thresholds"]["min_cells_per_sample"] # Min number of cells per sample when using pseudobulked data

SAMPLE_GROUP_COLS = [ c["colname"] for c in config["sample_group_pairs"] ]
SAMPLE_GROUP_LHSS = [ c["lhs"] for c in config["sample_group_pairs"] ]
SAMPLE_GROUP_RHSS = [ c["rhs"] for c in config["sample_group_pairs"] ]

UNIQUE_SAMPLE_GROUP_COLS = sorted(set(SAMPLE_GROUP_COLS))


# Rules
rule all:
  input:
    CLEANED_H5AD_PATH,
    join(INTERMEDIATE_DIR, "combined.subclass_l1.specimen.sum.pdata.h5ad")

# TODO: rules for doing diff exp tests, for either .adata.h5ad or .pdata.h5ad files
# This rule will produce the outputs required by the "all" rule.


# rule normalize_basic:
#   input:
#     join_zdone(ZARR_PATH, "uns", "comparison_metadata")
#   output:
#     join_zdone(ZARR_PATH, "uns", "comparison_metadata.normalize_basic")
#   resources:
#     slurm_partition="short",
#     runtime=60*2, # 2 hours
#     mem_mb=240_000, # 120 GB
#     cpus_per_task=2
#   shell:
#     """
#     compasce \
#         --zarr-path {ZARR_PATH} \
#         --function-name "normalize_basic"
#     """

# Reference: https://www.sc-best-practices.org/conditions/differential_gene_expression.html#pseudobulk

rule combine_splits:
  input:
    lambda w: expand(
      join(INTERMEDIATE_DIR, "{{cell_type_col}}.{cell_type_name_norm}.{{sample_id_col}}.{sample_id}.{{agg_func}}.agg.adata.h5ad"),
      cell_type_name_norm=[normalize_cell_type(ct) for ct in config["cell_types"][w.cell_type_col]],
      sample_id=SPECIMEN_IDS,
    )
  output:
    join(INTERMEDIATE_DIR, "combined.{cell_type_col}.{sample_id_col}.{agg_func}.pdata.h5ad")
  script:
    join(SCRIPTS_DIR, "combine_splits.py")


rule aggregate_split:
  input:
    join(INTERMEDIATE_DIR, "{cell_type_col}.{cell_type_name_norm}.{sample_id_col}.{sample_id}.filtered.adata.h5ad")
  output:
    join(INTERMEDIATE_DIR, "{cell_type_col}.{cell_type_name_norm}.{sample_id_col}.{sample_id}.{agg_func}.agg.adata.h5ad")
  params:
    cell_type_name_orig=lambda w: unnormalize_cell_type(w.cell_type_name_norm)
  shell:
    """
    python scripts/agg_adata.py \
        --input-h5ad {input} \
        --output-h5ad {output} \
        --cell-type-col {wildcards.cell_type_col} \
        --cell-type-name {params.cell_type_name_orig} \
        --sample-id-col {wildcards.sample_id_col} \
        --sample-id {wildcards.sample_id} \
        --agg-func {wildcards.agg_func}
    """

# Take a map-reduce approach to pseudobulking.
# In parallel, subset to each cell type and sample group as needed for the pseudobulk.
# Do not yet aggregate, but this can be trivially done in follow-up steps that run in parallel.
# In the individual split files, we could potentially save as dense feasibly as well.
rule split_for_pseudobulk_by_cell_type_and_specimen_id:
  input:
    CLEANED_H5AD_PATH
  output:
    join(INTERMEDIATE_DIR, "{cell_type_col}.{cell_type_name_norm}.{sample_id_col}.{sample_id}.filtered.adata.h5ad")
  params:
    cell_type_name_orig=lambda w: unnormalize_cell_type(w.cell_type_name_norm),
  shell:
    """
    python scripts/split_adata.py \
        --input-h5ad {CLEANED_H5AD_PATH} \
        --output-h5ad {output} \
        --cell-type-col {wildcards.cell_type_col} \
        --cell-type-name {params.cell_type_name_orig} \
        --sample-id-col {wildcards.sample_id_col} \
        --sample-id {wildcards.sample_id} \
    """

rule convert_to_zarr:
  input:
    h5ad=join(RAW_DIR, "kpmp-aug-2025", "SingleNucleus_KPMP_Explorer_05182025.h5ad"),
    clinical=join(RAW_DIR, "kpmp-aug-2025", "20250606_OpenAccessClinicalData.csv"),
    deg_dir=join(RAW_DIR, "kpmp-aug-2025")
  output:
    CLEANED_H5AD_PATH
  resources:
    slurm_partition="short",
    runtime=60*2, # 2 hours
    mem_mb=240_000, # 240 GB
    cpus_per_task=2
  shell:
    """
    python scripts/run_comparisons_kpmp_2025.py \
        --input-h5ad {input.h5ad} \
        --input-csv {input.clinical} \
        --input-deg-dir {input.deg_dir} \
        --output {CLEANED_H5AD_PATH} \
        --stop-early
    """

# No download rule:
# - Download raw data from Globus and put in RAW_DIR/kpmp-aug-2025
# - Download clinical metadata from  KPMP Atlas Repository
#   https://atlas.kpmp.org/repository/?size=n_20_n&filters%5B0%5D%5Bfield%5D=data_type&filters%5B0%5D%5Bvalues%5D%5B0%5D=Clinical%20Study%20Data&filters%5B0%5D%5Btype%5D=any
