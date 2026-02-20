# compasce-degs

Differential expression tests for KPMP single-nucleus RNA-seq data.

We will output results to tabular (CSV) format.
We can also disregard the need to convert the expression matrix to dense format for now, and can keep it sparse to save on memory during DEG computations.

We will work with H5AD files before producing tabular outputs.

Note that the dataset is over 20 GB, exceeding my laptop RAM of 16 GB, so we cannot forget about accounting for the out-of-memory problem at the computationally expensive steps that process the expression matrix, like computing pseudobulks.

We want to use Dask to enable computation on a laptop for individual computations, and we want to use Snakemake for parallelization.

For now, we will limit comptations to subclass L1 cell types (with potential to expand in the future).

We will use cell count and sample count thresholds to determine whether certain comparisons are feasible to run:
- at least 25 cells per participant
- at least 3 participants

We also want to record how many participants and cells-per-participant contributed to each comparison result/diff-exp test.

For all-samples cases (not split by sample groups), we can do both:
- pseudobulking via DESeq2/pyDeseq2
- findAllMarkers (non-pseudobulked) via Seurat

and we can allow users to toggle between these results in the UI.

When using the adjudicated sample categorizations, we want to aggregate the AKI subcategories of ATI and AIN (roll them up into a parent "AKI" category), since there are currently not enough samples for comparisons at the subcategory level to make sense.
