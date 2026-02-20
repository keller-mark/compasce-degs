include: "./common.smk"
configfile: "./scrnaseq_kpmp.yaml"

RAW_H5AD_PATH = join(RAW_DIR, "kpmp-aug-2025", "SingleNucleus_KPMP_Explorer_05182025.h5ad")
RAW_SAMPLES_PATH = join(RAW_DIR, "kpmp-aug-2025", "20250606_OpenAccessClinicalData.csv")

# Intermediate output paths
NORMALIZED_H5AD_PATH = join(INTERMEDIATE_DIR, "normalized.h5ad")

SAMPLE_GROUP_COLS = [ c["colname"] for c in config["sample_group_pairs"] ]
SAMPLE_GROUP_LHSS = [ c["lhs"] for c in config["sample_group_pairs"] ]
SAMPLE_GROUP_RHSS = [ c["rhs"] for c in config["sample_group_pairs"] ]

# Rules
rule all:
  input:
    # Expansions for L1 cell types, cell type vs rest, not pseudobulked
    expand(
      join(PROCESSED_DIR, "cell_type_vs_rest", "not_pseudobulked", "{cell_type_col}", "{cell_type}.csv"),
      cell_type_col=["subclass_l1"],
      cell_type=config["cell_types"]["subclass_l1"],
    ),
    # Expansions for L1 cell types, cell type vs rest, pseudobulked
    expand(
      join(PROCESSED_DIR, "cell_type_vs_rest", "pseudobulked", "{cell_type_col}", "{cell_type}.csv"),
      cell_type_col=["subclass_l1"],
      cell_type=config["cell_types"]["subclass_l1"],
    ),
    # Expansions for L1 cell types, pairwise comparisons betwen sample groups (group_col, lhs_group, rhs_group), pseudobulked
    expand(
      expand(
        join(PROCESSED_DIR, "within_cell_type", "pseudobulked", "{sample_col}", "{lhs_group}", "{rhs_group}", "{{cell_type_col}}", "{{cell_type}}.csv"),
        zip, # We use zip here to avoid combinatorial expansion, since we want to match the (sample_col, lhs_group, rhs_group) together as tuples.
        sample_col=SAMPLE_GROUP_COLS,
        lhs_group=SAMPLE_GROUP_LHSS,
        rhs_group=SAMPLE_GROUP_RHSS,
      ),
      cell_type_col="subclass_l1",
      cell_type=config["cell_types"]["subclass_l1"],
    )



# TODO: rules for doing diff exp tests, for either .adata.h5ad or .pdata.h5ad files
# This rule will produce the outputs required by the "all" rule.



# TODO: rules for computing pseudobulks, outputing as .pdata.h5ad files



rule normalize_basic:
  input:
    join_zdone(ZARR_PATH, "uns", "comparison_metadata")
  output:
    join_zdone(ZARR_PATH, "uns", "comparison_metadata.normalize_basic")
  resources:
    slurm_partition="short",
    runtime=60*2, # 2 hours
    mem_mb=240_000, # 120 GB
    cpus_per_task=2
  shell:
    """
    compasce \
        --zarr-path {ZARR_PATH} \
        --function-name "normalize_basic"
    """

rule convert_to_zarr:
  input:
    h5ad=join(RAW_DIR, "kpmp-aug-2025", "SingleNucleus_KPMP_Explorer_05182025.h5ad"),
    clinical=join(RAW_DIR, "kpmp-aug-2025", "20250606_OpenAccessClinicalData.csv"),
    deg_dir=join(RAW_DIR, "kpmp-aug-2025")
  output:
    join_zdone(ZARR_PATH, "uns", "comparison_metadata")
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
        --output {ZARR_PATH} \
        --stop-early
    """

# No download rule:
# - Download raw data from Globus and put in RAW_DIR/kpmp-aug-2025
# - Download clinical metadata from  KPMP Atlas Repository
#   https://atlas.kpmp.org/repository/?size=n_20_n&filters%5B0%5D%5Bfield%5D=data_type&filters%5B0%5D%5Bvalues%5D%5B0%5D=Clinical%20Study%20Data&filters%5B0%5D%5Btype%5D=any
