include: "./common.smk"
configfile: "./config.yaml"

RAW_H5AD_PATH = join(RAW_DIR, "kpmp-aug-2025", "SingleNucleus_KPMP_Explorer_05182025.h5ad")
RAW_SAMPLES_PATH = join(RAW_DIR, "kpmp-aug-2025", "20250606_OpenAccessClinicalData.csv")

# Intermediate output paths
NORMALIZED_H5AD_PATH = join(INTERMEDIATE_DIR, "normalized.h5ad")

# Build (sample_col, lhs_group, rhs_group) tuples from config for expand()
SAMPLE_GROUP_PAIRS = config["sample_group_pairs"]
SAMPLE_COLS   = [p["colname"] for p in SAMPLE_GROUP_PAIRS]
LHS_GROUPS    = [p["lhs"]     for p in SAMPLE_GROUP_PAIRS]
RHS_GROUPS    = [p["rhs"]     for p in SAMPLE_GROUP_PAIRS]

# Sanitize group names for use in file paths (spaces -> underscores, slashes -> _)
def sanitize(s):
    return s.replace(" ", "_").replace("/", "_")

SAMPLE_COLS_SAN  = [sanitize(c) for c in SAMPLE_COLS]
LHS_GROUPS_SAN   = [sanitize(g) for g in LHS_GROUPS]
RHS_GROUPS_SAN   = [sanitize(g) for g in RHS_GROUPS]

# Pseudobulk intermediate paths per (sample_col, cell_type_col, cell_type)
def pseudobulk_path(sample_col, cell_type_col, cell_type):
    return join(
        INTERMEDIATE_DIR, "pseudobulk",
        sanitize(sample_col), sanitize(cell_type_col), sanitize(cell_type) + ".h5ad"
    )


# Rules
rule all:
  input:
    # Cell type vs rest, not pseudobulked
    expand(
      join(PROCESSED_DIR, "cell_type_vs_rest", "not_pseudobulked", "{cell_type_col}", "{cell_type}.csv"),
      cell_type_col="subclass_l1",
      cell_type=config["cell_types"]["subclass_l1"],
    ),
    # Cell type vs rest, pseudobulked
    expand(
      join(PROCESSED_DIR, "cell_type_vs_rest", "pseudobulked", "{cell_type_col}", "{cell_type}.csv"),
      cell_type_col="subclass_l1",
      cell_type=config["cell_types"]["subclass_l1"],
    ),
    # Within cell type, sample group pairwise comparisons, pseudobulked
    # We zip (sample_col, lhs, rhs) to avoid the combinatorial explosion of expand()
    [
      expand(
        join(
            PROCESSED_DIR, "within_cell_type", "pseudobulked",
            SAMPLE_COLS_SAN[i], LHS_GROUPS_SAN[i], RHS_GROUPS_SAN[i],
            "{cell_type_col}", "{cell_type}.csv"
        ),
        cell_type_col="subclass_l1",
        cell_type=config["cell_types"]["subclass_l1"],
      )
      for i in range(len(SAMPLE_GROUP_PAIRS))
    ],


# ---------------------------------------------------------------------------
# Normalization
# ---------------------------------------------------------------------------

rule normalize_basic:
  input:
    h5ad=RAW_H5AD_PATH,
    csv=RAW_SAMPLES_PATH,
  output:
    NORMALIZED_H5AD_PATH
  shell:
    """
    python scripts/normalize.py \
      --input-h5ad {input.h5ad} \
      --input-csv  {input.csv} \
      --output     {output}
    """


# ---------------------------------------------------------------------------
# Pseudobulking  (one output per sample_col x cell_type_col x cell_type)
# We generate one rule per distinct sample_col to keep wildcards tractable.
# ---------------------------------------------------------------------------

rule pseudobulk:
  input:
    NORMALIZED_H5AD_PATH
  output:
    join(INTERMEDIATE_DIR, "pseudobulk", "{sample_col}", "{cell_type_col}", "{cell_type}.h5ad")
  wildcard_constraints:
    cell_type="[^/]+"
  shell:
    """
    python scripts/pseudobulk.py \
      --input          {input} \
      --output         {output} \
      --sample-col     {wildcards.sample_col} \
      --cell-type-col  {wildcards.cell_type_col} \
      --cell-type      "{wildcards.cell_type}"
    """


# ---------------------------------------------------------------------------
# DEG: cell type vs rest  (not pseudobulked, Seurat FindMarkers / Wilcoxon)
# ---------------------------------------------------------------------------

rule deg_cell_type_vs_rest_not_pseudobulked:
  input:
    NORMALIZED_H5AD_PATH
  output:
    join(PROCESSED_DIR, "cell_type_vs_rest", "not_pseudobulked", "{cell_type_col}", "{cell_type}.csv")
  wildcard_constraints:
    cell_type="[^/]+"
  shell:
    """
    Rscript scripts/deg_cell_type_vs_rest_nonpb.R \
      --input         {input} \
      --output        {output} \
      --cell-type-col {wildcards.cell_type_col} \
      --cell-type     "{wildcards.cell_type}"
    """


# ---------------------------------------------------------------------------
# DEG: cell type vs rest  (pseudobulked, pydeseq2)
# ---------------------------------------------------------------------------

rule deg_cell_type_vs_rest_pseudobulked:
  input:
    pb=join(INTERMEDIATE_DIR, "pseudobulk", "all_samples", "{cell_type_col}", "{cell_type}.h5ad"),
    norm=NORMALIZED_H5AD_PATH,
  output:
    join(PROCESSED_DIR, "cell_type_vs_rest", "pseudobulked", "{cell_type_col}", "{cell_type}.csv")
  wildcard_constraints:
    cell_type="[^/]+"
  shell:
    """
    python scripts/deg_cell_type_vs_rest_pb.py \
      --input           {input.pb} \
      --normalized-h5ad {input.norm} \
      --output          {output} \
      --cell-type-col   {wildcards.cell_type_col} \
      --cell-type       "{wildcards.cell_type}"
    """


# ---------------------------------------------------------------------------
# DEG: within cell type, sample group pairs (pseudobulked, pydeseq2)
# ---------------------------------------------------------------------------

rule deg_within_cell_type_pseudobulked:
  input:
    lambda wc: pseudobulk_path(wc.sample_col, wc.cell_type_col, wc.cell_type)
  output:
    join(PROCESSED_DIR, "within_cell_type", "pseudobulked",
         "{sample_col}", "{lhs_group}", "{rhs_group}", "{cell_type_col}", "{cell_type}.csv")
  wildcard_constraints:
    cell_type="[^/]+"
  shell:
    """
    python scripts/deg_within_cell_type_pb.py \
      --input         {input} \
      --output        {output} \
      --sample-col    {wildcards.sample_col} \
      --lhs-group     "{wildcards.lhs_group}" \
      --rhs-group     "{wildcards.rhs_group}" \
      --cell-type-col {wildcards.cell_type_col} \
      --cell-type     "{wildcards.cell_type}"
    """
