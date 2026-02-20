This is a bioinformatics repository, aiming to perform basic analysis of single-nucleus RNA-seq data:
- normalization
- pseudobulking
- multiple types of differential expression tests

Store raw and intermediate data using AnnData H5AD files.
Store outputs of differential expression tests using CSV files.

We are using the following tools:
- Snakemake for workflow management
- python and R
- scanpy (in python) and seurat (in R)
- pydeseq2 (in python) and deseq2 (in R)
- anndata (in python) and anndataR (for R interoperability of h5ad files)

Our dataset is very large (over 600,000 cells and 30,000 genes) and the full expression matrix does not fit in memory/RAM, requiring usage of solutions such as Dask and usage of Scipy sparse arrays.
