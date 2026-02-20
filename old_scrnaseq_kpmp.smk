include: "./common.smk"
configfile: "./scrnaseq_kpmp.yaml"

ZARR_PATH = join(PROCESSED_DIR, "kpmp-aug-2025.adata.zarr")

# Rules
rule all:
  input:
    join_zdone(ZARR_PATH, "uns", "comparison_metadata.merged")

rule merge_metadata:
  input:
    join_zdone(ZARR_PATH, "uns", "comparison_metadata"),
    join_zdone(ZARR_PATH, "uns", "comparison_metadata.normalize_basic"),
    join_zdone(ZARR_PATH, "uns", "comparison_metadata.normalize_pearson_residuals"),
    join_zdone(ZARR_PATH, "uns", "comparison_metadata.densmap"),
    #join_zdone(ZARR_PATH, "uns", "comparison_metadata.compute_diffexp"),
    join_zdone(ZARR_PATH, "uns", "comparison_metadata.compute_diffexp_pydeseq2"),
    join_zdone(ZARR_PATH, "uns", "comparison_metadata.compute_diffabundance"),
    #join_zdone(ZARR_PATH, "uns", "comparison_metadata.compute_lemur")
  output:
    join_zdone(ZARR_PATH, "uns", "comparison_metadata.merged")
  resources:
    slurm_partition="short",
    runtime=60, # 1 hour
    mem_mb=4_000, # 4 GB
    cpus_per_task=2
  shell:
    """
    python scripts/merge_metadata.py \
        --zarr-path {ZARR_PATH}
    """

rule compute_lemur:
  input:
    join_zdone(ZARR_PATH, "uns", "comparison_metadata.normalize_basic")
  output:
    join_zdone(ZARR_PATH, "uns", "comparison_metadata.compute_lemur")
  resources:
    slurm_partition="medium",
    runtime=60*24*4, # 4 days
    mem_mb=240_000, # 240 GB
    cpus_per_task=4
  shell:
    """
    compasce \
        --zarr-path {ZARR_PATH} \
        --function-name "compute_lemur" \
        --mem-limit "24GB" \
        --n-workers 4 \
        --threads-per-worker 2
    """

rule compute_diffexp:
  input:
    metadata_path=join_zdone(ZARR_PATH, "uns", "comparison_metadata.normalize_basic"),
    deg_dir=join(RAW_DIR, "kpmp-aug-2025")
  output:
    join_zdone(ZARR_PATH, "uns", "comparison_metadata.compute_diffexp")
  resources:
    slurm_partition="medium",
    runtime=60*24*5, # 5 days
    mem_mb=240_000, # 240 GB
    cpus_per_task=4
  shell:
    """
    compasce \
        --zarr-path {ZARR_PATH} \
        --input-deg-dir {input.deg_dir} \
        --function-name "compute_diffexp"
    """

rule compute_diffexp_pydeseq2:
  input:
    metadata_path=join_zdone(ZARR_PATH, "uns", "comparison_metadata.normalize_basic"),
    deg_dir=join(RAW_DIR, "kpmp-aug-2025")
  output:
    join_zdone(ZARR_PATH, "uns", "comparison_metadata.compute_diffexp_pydeseq2")
  resources:
    slurm_partition="medium",
    runtime=60*24*5, # 5 days
    mem_mb=240_000, # 240 GB
    cpus_per_task=4
  shell:
    """
    compasce \
        --zarr-path {ZARR_PATH} \
        --input-deg-dir {input.deg_dir} \
        --function-name "compute_diffexp_pydeseq2"
    """

rule compute_diffabundance:
  input:
    join_zdone(ZARR_PATH, "uns", "comparison_metadata.normalize_basic"),
    join_zdone(ZARR_PATH, "uns", "comparison_metadata.compute_diffexp_pydeseq2") # TEMP
  output:
    join_zdone(ZARR_PATH, "uns", "comparison_metadata.compute_diffabundance")
  resources:
    slurm_partition="medium",
    runtime=60*24*2, # 2 days
    mem_mb=120_000, # 120 GB
    cpus_per_task=2
  shell:
    """
    compasce \
        --zarr-path {ZARR_PATH} \
        --function-name "compute_diffabundance"
    """

rule densmap:
  input:
    join_zdone(ZARR_PATH, "uns", "comparison_metadata.normalize_basic"),
    join_zdone(ZARR_PATH, "uns", "comparison_metadata.compute_diffexp_pydeseq2") # TEMP
  output:
    join_zdone(ZARR_PATH, "uns", "comparison_metadata.densmap")
  resources:
    slurm_partition="short",
    runtime=60*11, # 11 hours
    mem_mb=240_000, # 240 GB
    cpus_per_task=4
  shell:
    """
    compasce \
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
    compasce \
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
