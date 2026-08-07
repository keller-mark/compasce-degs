include: "./common.smk"
configfile: "./scrnaseq_kpmp.yaml"

DEBUG_MODE = False

RAW_H5AD_PATH = join(RAW_DIR, "kpmp-aug-2025", "SingleNucleus_KPMP_Explorer_05182025.h5ad")
RAW_SAMPLES_PATH = join(RAW_DIR, "kpmp-aug-2025", "20250606_OpenAccessClinicalData.csv")

# Intermediate output paths
CLEANED_H5AD_PATH = join(INTERMEDIATE_DIR, "cleaned.h5ad")
ZARR_PATH = join(PROCESSED_DIR, "kpmp-apr-2026.adata.zarr")
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
  
def normalize_identifier(cell_type):
    # Normalize cell type names to be filesystem-friendly, since we'll be using them in file paths.
    return (
      cell_type
        .replace("/", "__slash__")
        .replace(" ", "__space__")
        .replace("+", "__plus__")
        .replace("\xEF", "__i__")
    )

def unnormalize_identifier(cell_type):
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
NUM_COUNTS_PER_GENE_THRESHOLD = config["thresholds"]["min_counts_per_gene"] # Min number of counts per gene (post-pseudobulk-aggregation) for a gene to be considered expressed in the pseudobulk sample (independently computed per cell type).
FRAC_SAMPLES_PER_GENE_THRESHOLD = config["thresholds"]["min_frac_samples_per_gene"] # Min percentage of samples expressing gene (post-pseudobulk-aggregation) for a gene to be included in comparisons.

SAMPLE_GROUP_COLS = [ c["colname"] for c in config["sample_group_pairs"] ]
SAMPLE_GROUP_LHSS = [ c["lhs"] for c in config["sample_group_pairs"] ]
SAMPLE_GROUP_RHSS = [ c["rhs"] for c in config["sample_group_pairs"] ]

UNIQUE_SAMPLE_GROUP_COLS = sorted(set(SAMPLE_GROUP_COLS))


if DEBUG_MODE:
    L1_CELL_TYPES = L1_CELL_TYPES[:5]
    L2_CELL_TYPES = L2_CELL_TYPES[:5]
    L3_CELL_TYPES = L3_CELL_TYPES[:5]
    SPECIMEN_IDS = SPECIMEN_IDS[:5]


# Rules
rule all:
  input:
    CLEANED_H5AD_PATH,
    join(INTERMEDIATE_DIR, "combined.subclass_l1.specimen.sum.pdata.h5ad"),
    
    join_zdone(ZARR_PATH, "uns", "comparison_metadata.normalize_basic"),
    join_zdone(ZARR_PATH, "uns", "comparison_metadata.normalize_pearson_residuals"),
    join_zdone(ZARR_PATH, "uns", "comparison_metadata.densmap"),
    join_zdone(ZARR_PATH, "uns", "comparison_metadata.compute_diffabundance"),

    # Subclass L1
    # L1: Cell type vs rest
    expand(
      #join(INTERMEDIATE_DIR, "pydeseq.{cell_type_col}.{cell_type_name_norm}.{sample_id_col}.{agg_func}.csv"),
      
      join_zdone(ZARR_PATH, "uns", "comparison_metadata.pydeseq.{cell_type_col}.{cell_type_name_norm}.{sample_id_col}.{agg_func}"),
      cell_type_col=["subclass_l1"],
      cell_type_name_norm=[normalize_identifier(ct) for ct in L1_CELL_TYPES],
      sample_id_col=[SPECIMEN_ID_COL],
      agg_func=["sum"]
    ),
    # L1: Within cell type, case vs control
    expand(
      [
        #join(INTERMEDIATE_DIR, f"pydeseq_within_celltype.{{cell_type_col}}.{{cell_type_name_norm}}.{{sample_id_col}}.{c['colname']}.{normalize_identifier(c['lhs'])}.{normalize_identifier(c['rhs'])}.{{agg_func}}.csv")
        
        #join_zdone(ZARR_PATH, "uns", f"comparison_metadata.pydeseq_within_celltype.{{cell_type_col}}.{{cell_type_name_norm}}.{{sample_id_col}}.{c['colname']}.{normalize_identifier(c['lhs'])}.{normalize_identifier(c['rhs'])}.{{agg_func}}")
        join_zdone(ZARR_PATH, "uns", "comparison_metadata.pydeseq_within_celltype.{cell_type_col}.{cell_type_name_norm}.{sample_id_col}." + c['colname'] + "." + normalize_identifier(c['lhs']) + "." + normalize_identifier(c['rhs']) + ".{agg_func}")
        for c in config["sample_group_pairs"]
      ],
      cell_type_col=["subclass_l1"],
      cell_type_name_norm=[normalize_identifier(ct) for ct in L1_CELL_TYPES],
      sample_id_col=[SPECIMEN_ID_COL],
      agg_func=["sum"]
    )
    # TODO: uncomment
    # ,
    # # Subclass L2
    # # L2: Cell type vs rest
    # expand(
    #   join(INTERMEDIATE_DIR, "pydeseq.{cell_type_col}.{cell_type_name_norm}.{sample_id_col}.{agg_func}.csv"),
    #   cell_type_col=["subclass_l2"],
    #   cell_type_name_norm=[normalize_identifier(ct) for ct in L2_CELL_TYPES],
    #   sample_id_col=[SPECIMEN_ID_COL],
    #   agg_func=["sum"]
    # ),
    # # L2: Within cell type, case vs control
    # expand(
    #   [
    #     join(INTERMEDIATE_DIR, f"pydeseq_within_celltype.{{cell_type_col}}.{{cell_type_name_norm}}.{{sample_id_col}}.{c['colname']}.{normalize_identifier(c['lhs'])}.{normalize_identifier(c['rhs'])}.{{agg_func}}.csv")
    #     for c in config["sample_group_pairs"]
    #   ],
    #   cell_type_col=["subclass_l2"],
    #   cell_type_name_norm=[normalize_identifier(ct) for ct in L2_CELL_TYPES],
    #   sample_id_col=[SPECIMEN_ID_COL],
    #   agg_func=["sum"]
    # )


# Convert H5AD to Zarr, normalize, run densmap.
# Then, insert the pydeseq output dataframes into the Zarr store, along with the comparison metadata.



# RULES TO INSERT PYDESEQ RESULTS INTO ZARR STORE WITH METADATA

# TODO: Insert these dataframes into the smaller PSEUDOBULKED anndata object, rather than the single-cell-resolution anndata object?

rule insert_within_celltype_case_vs_control_degs:
  input:
    ladata=join_zdone(ZARR_PATH, "uns", "comparison_metadata"),
    deg_results=join(INTERMEDIATE_DIR, "pydeseq_within_celltype.{cell_type_col}.{cell_type_name_norm}.{sample_id_col}.{sample_group_col}.{sample_group_lhs_norm}.{sample_group_rhs_norm}.{agg_func}.csv")
  output:
    join_zdone(ZARR_PATH, "uns", "comparison_metadata.pydeseq_within_celltype.{cell_type_col}.{cell_type_name_norm}.{sample_id_col}.{sample_group_col}.{sample_group_lhs_norm}.{sample_group_rhs_norm}.{agg_func}")
  params:
    cell_type_name_orig=lambda w: unnormalize_identifier(w.cell_type_name_norm),
    sample_group_lhs_orig=lambda w: unnormalize_identifier(w.sample_group_lhs_norm),
    sample_group_rhs_orig=lambda w: unnormalize_identifier(w.sample_group_rhs_norm)
  resources:
    slurm_partition="short",
    runtime=60, # half hour
    mem_mb=16_000, # 16 GB
    cpus_per_task=2
  shell:
    """
    compasce2 \
        --zarr-path {ZARR_PATH} \
    insert_within_celltype_case_vs_control_degs \
        --csv-path {input.deg_results} \
        --cell-type-col {wildcards.cell_type_col} \
        --cell-type-name "{params.cell_type_name_orig}" \
        --sample-id-col {wildcards.sample_id_col} \
        --sample-group-col "{wildcards.sample_group_col}" \
        --sample-group-lhs "{params.sample_group_lhs_orig}" \
        --sample-group-rhs "{params.sample_group_rhs_orig}" \
        --agg-func {wildcards.agg_func} \
        --out-path {output}
    """

rule insert_celltype_vs_rest_degs:
  input:
    ladata=join_zdone(ZARR_PATH, "uns", "comparison_metadata"),
    deg_results=join(INTERMEDIATE_DIR, "pydeseq.{cell_type_col}.{cell_type_name_norm}.{sample_id_col}.{agg_func}.csv")
  output:
    join_zdone(ZARR_PATH, "uns", "comparison_metadata.pydeseq.{cell_type_col}.{cell_type_name_norm}.{sample_id_col}.{agg_func}")
  params:
    cell_type_name_orig=lambda w: unnormalize_identifier(w.cell_type_name_norm)
  resources:
    slurm_partition="short",
    runtime=60, # half hour
    mem_mb=16_000, # 16 GB
    cpus_per_task=2
  shell:
    """
    compasce2 \
        --zarr-path {ZARR_PATH} \
    insert_celltype_vs_rest_degs \
        --csv-path {input.deg_results} \
        --cell-type-col {wildcards.cell_type_col} \
        --cell-type-name "{params.cell_type_name_orig}" \
        --sample-id-col {wildcards.sample_id_col} \
        --agg-func {wildcards.agg_func} \
        --out-path {output}
    """


# Begin rules copied from https://github.com/keller-mark/compasce/blob/keller-mark/kpmp-nov-2025/scrnaseq_kpmp.smk
rule compute_diffabundance:
  input:
    join_zdone(ZARR_PATH, "uns", "comparison_metadata.normalize_basic")
  output:
    join_zdone(ZARR_PATH, "uns", "comparison_metadata.compute_diffabundance")
  resources:
    slurm_partition="medium",
    runtime=60*24*2, # 2 days
    mem_mb=120_000, # 120 GB
    cpus_per_task=2
  shell:
    """
    compasce2 \
        --zarr-path {ZARR_PATH} \
        --function-name "compute_diffabundance"
    """

rule densmap:
  input:
    join_zdone(ZARR_PATH, "uns", "comparison_metadata.normalize_basic")
  output:
    join_zdone(ZARR_PATH, "uns", "comparison_metadata.densmap")
  resources:
    slurm_partition="short",
    runtime=60*11, # 11 hours
    mem_mb=240_000, # 240 GB
    cpus_per_task=4
  shell:
    """
    compasce2 \
        --zarr-path {ZARR_PATH} \
        --function-name "densmap"
    """

rule normalize_pearson_residuals:
  input:
    join_zdone(ZARR_PATH, "uns", "comparison_metadata")
  output:
    join_zdone(ZARR_PATH, "uns", "comparison_metadata.normalize_pearson_residuals")
  resources:
    slurm_partition="short",
    runtime=60*2, # 2 hours
    mem_mb=120_000, # 120 GB
    cpus_per_task=2
  shell:
    """
    compasce2 \
        --zarr-path {ZARR_PATH} \
        --function-name "normalize_pearson_residuals"
    """

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
    compasce2 \
        --zarr-path {ZARR_PATH} \
        --function-name "normalize_basic"
    """

# End rules copied from https://github.com/keller-mark/compasce/blob/keller-mark/kpmp-nov-2025/scrnaseq_kpmp.smk




# Reference: https://www.sc-best-practices.org/conditions/differential_gene_expression.html#pseudobulk

rule pydeseq_within_celltype_case_vs_control:
  input:
    join(INTERMEDIATE_DIR, "combined.{cell_type_col}.{sample_id_col}.{agg_func}.pdata.h5ad")
  output:
    de_df=join(INTERMEDIATE_DIR, "pydeseq_within_celltype.{cell_type_col}.{cell_type_name_norm}.{sample_id_col}.{sample_group_col}.{sample_group_lhs_norm}.{sample_group_rhs_norm}.{agg_func}.csv"),
    obs_filtering_df=join(INTERMEDIATE_DIR, "pydeseq_within_celltype_obs_filtering.{cell_type_col}.{cell_type_name_norm}.{sample_id_col}.{sample_group_col}.{sample_group_lhs_norm}.{sample_group_rhs_norm}.{agg_func}.csv"),
    var_filtering_df=join(INTERMEDIATE_DIR, "pydeseq_within_celltype_var_filtering.{cell_type_col}.{cell_type_name_norm}.{sample_id_col}.{sample_group_col}.{sample_group_lhs_norm}.{sample_group_rhs_norm}.{agg_func}.csv")
  params:
    cell_type_name_orig=lambda w: unnormalize_identifier(w.cell_type_name_norm),
    sample_group_lhs_orig=lambda w: unnormalize_identifier(w.sample_group_lhs_norm),
    sample_group_rhs_orig=lambda w: unnormalize_identifier(w.sample_group_rhs_norm)
  resources:
    slurm_partition="short",
    runtime=60, # 1 hour
    mem_mb=32_000, # 16 GB
    cpus_per_task=2
  shell:
    """
    python scripts/60_pydeseq_within_celltype_case_vs_control.py \
        --input-h5ad {input} \
        --output-de-csv {output.de_df} \
        --output-obs-filtering-csv {output.obs_filtering_df} \
        --output-var-filtering-csv {output.var_filtering_df} \
        --cell-type-col {wildcards.cell_type_col} \
        --cell-type-name "{params.cell_type_name_orig}" \
        --sample-id-col {wildcards.sample_id_col} \
        --sample-group-col "{wildcards.sample_group_col}" \
        --sample-group-lhs "{params.sample_group_lhs_orig}" \
        --sample-group-rhs "{params.sample_group_rhs_orig}" \
        --num-samples-threshold {NUM_SAMPLES_THRESHOLD} \
        --num-cells-per-sample-threshold {NUM_CELLS_PER_SAMPLE_THRESHOLD} \
        --num-counts-per-gene-threshold {NUM_COUNTS_PER_GENE_THRESHOLD} \
        --frac-samples-per-gene-threshold {FRAC_SAMPLES_PER_GENE_THRESHOLD}
    """

rule pydeseq_celltype_vs_rest:
  input:
    join(INTERMEDIATE_DIR, "combined.{cell_type_col}.{sample_id_col}.{agg_func}.pdata.h5ad")
  output:
    de_df=join(INTERMEDIATE_DIR, "pydeseq.{cell_type_col}.{cell_type_name_norm}.{sample_id_col}.{agg_func}.csv"),
    obs_filtering_df=join(INTERMEDIATE_DIR, "pydeseq_obs_filtering.{cell_type_col}.{cell_type_name_norm}.{sample_id_col}.{agg_func}.csv"),
    var_filtering_df=join(INTERMEDIATE_DIR, "pydeseq_var_filtering.{cell_type_col}.{cell_type_name_norm}.{sample_id_col}.{agg_func}.csv")
  params:
    cell_type_name_orig=lambda w: unnormalize_identifier(w.cell_type_name_norm)
  resources:
    slurm_partition="short",
    runtime=60, # 1 hour
    mem_mb=32_000, # 16 GB
    cpus_per_task=2
  shell:
    """
    python scripts/50_pydeseq_celltype_vs_rest.py \
        --input-h5ad {input} \
        --output-de-csv {output.de_df} \
        --output-obs-filtering-csv {output.obs_filtering_df} \
        --output-var-filtering-csv {output.var_filtering_df} \
        --sample-id-col {wildcards.sample_id_col} \
        --cell-type-col {wildcards.cell_type_col} \
        --cell-type-name "{params.cell_type_name_orig}" \
        --num-samples-threshold {NUM_SAMPLES_THRESHOLD} \
        --num-cells-per-sample-threshold {NUM_CELLS_PER_SAMPLE_THRESHOLD} \
        --num-counts-per-gene-threshold {NUM_COUNTS_PER_GENE_THRESHOLD} \
        --frac-samples-per-gene-threshold {FRAC_SAMPLES_PER_GENE_THRESHOLD}
    """

rule combine_splits:
  input:
    lambda w: expand(
      join(INTERMEDIATE_DIR, "{{cell_type_col}}.{cell_type_name_norm}.{{sample_id_col}}.{sample_id}.{{agg_func}}.agg.adata.h5ad"),
      cell_type_name_norm=[normalize_identifier(ct) for ct in config["cell_types"][w.cell_type_col]],
      sample_id=SPECIMEN_IDS,
    )
  output:
    protected(join(INTERMEDIATE_DIR, "combined.{cell_type_col}.{sample_id_col}.{agg_func}.pdata.h5ad"))
  resources:
    slurm_partition="short",
    runtime=30, # half hour
    mem_mb=32_000, # 16 GB
    cpus_per_task=2
  script:
    join(SCRIPTS_DIR, "30_combine_splits.py")


rule aggregate_split:
  input:
    join(INTERMEDIATE_DIR, "{cell_type_col}.{cell_type_name_norm}.{sample_id_col}.{sample_id}.filtered.adata.h5ad")
  output:
    temp(join(INTERMEDIATE_DIR, "{cell_type_col}.{cell_type_name_norm}.{sample_id_col}.{sample_id}.{agg_func}.agg.adata.h5ad"))
  params:
    cell_type_name_orig=lambda w: unnormalize_identifier(w.cell_type_name_norm)
  resources:
    slurm_partition="short",
    runtime=30, # half hour
    mem_mb=16_000, # 16 GB
    cpus_per_task=2
  shell:
    """
    python scripts/20_agg_adata.py \
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
# (In the individual split files, we could potentially save as dense feasibly as well.)
rule split_for_pseudobulk_by_cell_type_and_specimen_id:
  input:
    CLEANED_H5AD_PATH
  output:
    temp(join(INTERMEDIATE_DIR, "{cell_type_col}.{cell_type_name_norm}.{sample_id_col}.{sample_id}.filtered.adata.h5ad"))
  params:
    cell_type_name_orig=lambda w: unnormalize_identifier(w.cell_type_name_norm),
  resources:
    slurm_partition="short",
    runtime=30, # half hour
    mem_mb=16_000, # 16 GB
    cpus_per_task=2
  shell:
    """
    python scripts/10_split_adata.py \
        --input-h5ad {CLEANED_H5AD_PATH} \
        --output-h5ad {output} \
        --cell-type-col {wildcards.cell_type_col} \
        --cell-type-name {params.cell_type_name_orig} \
        --sample-id-col {wildcards.sample_id_col} \
        --sample-id {wildcards.sample_id} \
    """



rule clean_h5ad:
  input:
    h5ad=join(RAW_DIR, "kpmp-aug-2025", "SingleNucleus_KPMP_Explorer_05182025.h5ad"),
    clinical=join(RAW_DIR, "kpmp-aug-2025", "20250606_OpenAccessClinicalData.csv"),
    deg_dir=join(RAW_DIR, "kpmp-aug-2025")
  output:
    protected(CLEANED_H5AD_PATH),
    join_zdone(ZARR_PATH, "uns", "comparison_metadata")
  params:
    subset_line=("--subset" if DEBUG_MODE else "--no-subset")
  resources:
    slurm_partition="short",
    runtime=60*2, # 2 hours
    mem_mb=240_000, # 240 GB
    cpus_per_task=2
  shell:
    """
    python scripts/00_run_comparisons_kpmp_2025.py \
        --input-h5ad {input.h5ad} \
        --input-csv {input.clinical} \
        --input-deg-dir {input.deg_dir} \
        --output {CLEANED_H5AD_PATH} \
        --output-zarr {ZARR_PATH} \
        --stop-early {params.subset_line}
    """

# No download rule:
# - Download raw data from Globus and put in RAW_DIR/kpmp-aug-2025
# - Download clinical metadata from  KPMP Atlas Repository
#   https://atlas.kpmp.org/repository/?size=n_20_n&filters%5B0%5D%5Bfield%5D=data_type&filters%5B0%5D%5Bvalues%5D%5B0%5D=Clinical%20Study%20Data&filters%5B0%5D%5Btype%5D=any
