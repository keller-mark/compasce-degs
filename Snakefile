include: "./common.smk"
configfile: "./config.yaml"

RAW_H5AD_PATH = join(RAW_DIR, "kpmp-aug-2025", "SingleNucleus_KPMP_Explorer_05182025.h5ad")
RAW_SAMPLES_PATH = join(RAW_DIR, "kpmp-aug-2025", "20250606_OpenAccessClinicalData.csv")

# Intermediate output paths
NORMALIZED_H5AD_PATH = join(INTERMEDIATE_DIR, "normalized.h5ad")

# Rules
rule all:
  input:
    # Expansions for L1 cell types, cell type vs rest, not pseudobulked
    expand(
      join(PROCESSED_DIR, "cell_type_vs_rest", "not_pseudobulked", "{cell_type_col}", "{cell_type}.csv"),
      cell_type_col="subclass_l1",
      cell_type=config["cell_types"]["subclass_l1"],
    ),
    # Expansions for L1 cell types, cell type vs rest, pseudobulked
    expand(
      join(PROCESSED_DIR, "cell_type_vs_rest", "pseudobulked", "{cell_type_col}", "{cell_type}.csv"),
      cell_type_col="subclass_l1",
      cell_type=config["cell_types"]["subclass_l1"],
    ),
    # Expansions for L1 cell types, pairwise comparisons betwen sample groups (group_col, lhs_group, rhs_group), pseudobulked
    expand(
      join(PROCESSED_DIR, "within_cell_type", "pseudobulked", "{sample_col}", "{lhs_group}", "{rhs_group}", "{cell_type_col}", "{cell_type}.csv"),
      cell_type_col="subclass_l1",
      cell_type=config["cell_types"]["subclass_l1"],
      # TODO: fill in the rest of wildcards
    ),



# TODO: rules for doing diff exp tests, for either .adata.h5ad or .pdata.h5ad files
# This rule will produce the outputs required by the "all" rule.



# TODO: rules for computing pseudobulks, outputing as .pdata.h5ad files



rule normalize_basic:
  input:
    RAW_H5AD_PATH
  output:
    NORMALIZED_H5AD_PATH
  shell:
    """
    # TODO
    """

