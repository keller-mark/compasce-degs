# deg_cell_type_vs_rest_nonpb.R
#
# Differential expression: one cell type vs. all other cells (not pseudobulked).
#
# Uses Seurat's FindMarkers with the Wilcoxon rank-sum test.
# Input is the normalized H5AD file (produced by normalize.py).
# The log1p_norm layer is used as the expression matrix for the test.
#
# Output CSV columns
# ------------------
# gene            gene name / ID
# avg_log2FC      average log2 fold change (Seurat convention)
# pval            uncorrected p-value  (p_val in Seurat output)
# pval_adj        Bonferroni-corrected p-value  (p_val_adj)
# pct_fg          fraction of cells in foreground expressing the gene  (pct.1)
# pct_bg          fraction of cells in background expressing the gene  (pct.2)
# n_cells_fg      number of cells in the foreground (cell type of interest)
# n_cells_bg      number of cells in the background (all other cells)
# n_donors_fg     number of unique donors in the foreground
# n_donors_bg     number of unique donors in the background
# cell_type       the target cell type
# cell_type_col   the obs column used for grouping

suppressPackageStartupMessages({
  library(anndataR)
  library(Seurat)
  library(dplyr)
})

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(args) {
  out <- list()
  i <- 1
  while (i <= length(args)) {
    key <- sub("^--", "", args[i])
    key <- gsub("-", "_", key)  # --cell-type-col -> cell_type_col
    val <- args[i + 1]
    out[[key]] <- val
    i <- i + 2
  }
  out
}

opts <- parse_args(args)

required <- c("input", "output", "cell_type_col", "cell_type")
missing_args <- setdiff(required, names(opts))
if (length(missing_args) > 0) {
  stop(paste("Missing required arguments:", paste(paste0("--", gsub("_", "-", missing_args)), collapse = ", ")))
}

input_path    <- opts$input
output_path   <- opts$output
cell_type_col <- opts$cell_type_col
cell_type_san <- opts$cell_type   # may be sanitized (spaces/slashes replaced)

# ---------------------------------------------------------------------------
# Helper: recover original value from a sanitized wildcard
# ---------------------------------------------------------------------------
unsanitize <- function(sanitized, candidates) {
  normed <- gsub("/", "_", gsub(" ", "_", candidates))
  idx <- which(normed == sanitized)
  if (length(idx) > 0) candidates[idx[1]] else sanitized
}

# ---------------------------------------------------------------------------
# Load H5AD as AnnData first to resolve the cell type label
# ---------------------------------------------------------------------------
cat(sprintf("Loading: %s\n", input_path))
adata <- read_h5ad(input_path)

# Recover actual cell type label (handles VSM/P -> VSM_P etc.)
all_cell_types <- unique(adata$obs[[cell_type_col]])
cell_type <- unsanitize(cell_type_san, all_cell_types)
if (cell_type != cell_type_san) {
  cat(sprintf("  Resolved cell type '%s' -> '%s'\n", cell_type_san, cell_type))
}

donor_col <- "patient"

# Compute counts before converting to Seurat (obs is an R data.frame in anndataR)
obs <- adata$obs
fg_cells <- obs[[cell_type_col]] == cell_type
n_fg <- sum(fg_cells)
n_bg <- sum(!fg_cells)

if (n_fg == 0) {
  stop(sprintf("No cells found for cell type '%s' in column '%s'.", cell_type, cell_type_col))
}

n_donors_fg <- if (donor_col %in% colnames(obs)) {
  length(unique(obs[[donor_col]][fg_cells]))
} else NA_integer_

n_donors_bg <- if (donor_col %in% colnames(obs)) {
  length(unique(obs[[donor_col]][!fg_cells]))
} else NA_integer_

cat(sprintf("  Foreground: %d cells (%s donors), Background: %d cells (%s donors)\n",
            n_fg,
            ifelse(is.na(n_donors_fg), "NA", as.character(n_donors_fg)),
            n_bg,
            ifelse(is.na(n_donors_bg), "NA", as.character(n_donors_bg))))

# ---------------------------------------------------------------------------
# Convert to Seurat, mapping log1p_norm layer as X (-> data slot)
# ---------------------------------------------------------------------------
cat("Converting to Seurat object...\n")

# x_mapping = "log1p_norm" tells as_Seurat to use the log1p_norm layer as X,
# which populates the 'data' slot that FindMarkers (wilcox) reads from.
# We also bring counts along so Seurat is fully formed.
seurat_obj <- adata$as_Seurat(
  x_mapping      = "log1p_norm",
  layers_mapping = c("counts")
)

# ---------------------------------------------------------------------------
# Set identity classes and run FindMarkers
# ---------------------------------------------------------------------------
Idents(seurat_obj) <- seurat_obj@meta.data[[cell_type_col]]

cat(sprintf("Running FindMarkers (Wilcoxon) for '%s' vs rest...\n", cell_type))

markers <- FindMarkers(
  seurat_obj,
  ident.1         = cell_type,
  ident.2         = NULL,   # NULL = all other cells
  test.use        = "wilcox",
  logfc.threshold = 0,      # return all genes, filter later
  min.pct         = 0,
  only.pos        = FALSE,
  verbose         = FALSE
)

# ---------------------------------------------------------------------------
# Tidy output
# ---------------------------------------------------------------------------
result <- markers %>%
  tibble::rownames_to_column("gene") %>%
  rename(
    pval     = p_val,
    pval_adj = p_val_adj,
    pct_fg   = pct.1,
    pct_bg   = pct.2,
  ) %>%
  mutate(
    n_cells_fg    = as.integer(n_fg),
    n_cells_bg    = as.integer(n_bg),
    n_donors_fg   = as.integer(n_donors_fg),
    n_donors_bg   = as.integer(n_donors_bg),
    cell_type     = cell_type,
    cell_type_col = cell_type_col,
  ) %>%
  select(gene, avg_log2FC, pval, pval_adj, pct_fg, pct_bg,
         n_cells_fg, n_cells_bg, n_donors_fg, n_donors_bg,
         cell_type, cell_type_col)

# ---------------------------------------------------------------------------
# Write output
# ---------------------------------------------------------------------------
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
write.csv(result, output_path, row.names = FALSE)
cat(sprintf("Wrote %d rows to %s\n", nrow(result), output_path))
