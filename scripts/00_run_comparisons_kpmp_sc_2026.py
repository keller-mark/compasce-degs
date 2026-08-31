# Single-cell
from compasce_degs import run_all, create_dask_client, create_o2_dask_client
from anndata import read_h5ad
import numpy as np
import pandas as pd
import argparse
import h5py


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-h5ad", type=str, required=True, help = "Path to KPMP H5AD file from the Atlas, august 2026.")
    parser.add_argument("--input-csv", type=str, required=True, help = "Path to KPMP clinical data CSV file.")
    parser.add_argument("--output", type=str, required=True, help = "Path to output H5AD file")
    parser.add_argument("--output-zarr", type=str, required=True, help = "Path to output Zarr directory")
    parser.add_argument("--subset", action=argparse.BooleanOptionalAction, default=False)
    parser.add_argument("--mem-limit", type=str, default='16GB', required=False)
    parser.add_argument("--overwrite", action=argparse.BooleanOptionalAction, default=False)
    parser.add_argument("--stop-early", action=argparse.BooleanOptionalAction, default=False)
    args = parser.parse_args()



    def get_adata():
        adata = read_h5ad(args.input_h5ad)

        should_subset = args.subset
        if should_subset:
            print("SUBSETTING")
            # subset using random sample so that multiple sample groups are represented to enable comparison
            np.random.seed(1)
            obs_subset = np.random.choice(adata.obs.index.tolist(), size=20_000, replace=False).tolist()
            var_slice = slice(None)
            adata = adata[obs_subset, var_slice].copy()
        else:
            print("NOT SUBSETTING")

        # Join adata.obs with clinical data from CSV
        clinical_data = pd.read_csv(args.input_csv)

        adata.obs = adata.obs.merge(clinical_data, left_on="KPMP Participant ID", right_on="Participant ID", how="left")

        # We could have done a left join, but then we would have to filter out samples that do not have clinical data later.
        # We also do not want strings to be converted to NaN, as these cause Zarr writing errors like "TypeError: expected unicode string, found nan".

        # This effectively does an inner join. We cannot use how="inner", since this would only affect adata.obs, and not other anndata fields.
        has_clinical_data = ~adata.obs["Participant ID"].isna()
        adata = adata[has_clinical_data, :].copy()

        print(adata.obs.head())

        adata.obs["Primary Adjudicated Category"] = adata.obs["Primary Adjudicated Category"].fillna("NA")

        # Cleanup of Enrollment Category for healthy reference to exclude certain Sample Type values.
        def clean_enrollment_category(row):
            if row["Enrollment Category"] == "Healthy Reference":
                healthy_sample_type = row["Sample Type"]
                if healthy_sample_type in ["Tumor Nephrectomy", "Total Tumor Nephrectomy", "Partial Tumor Nephrectomy", "Deceased Donor Nephrectomy"]:
                    # Append the Sample Type, which will result in this row not matching the plain "Healthy Reference" class.
                    return f"{row['Enrollment Category']} - {row['Sample Type']}"
                # The other healthy sample types are: "Transplant Pre-perfusion Biopsy", "Intra-operative Biopsy" (HRT Percutaneous Nephrolithotomy Protocol i.e. kidney stone).
                # We want these types of healthy sample types to correspond to "Healthy Reference" (no suffix).
            return row["Enrollment Category"]

        # Cleanup of sample-level data
        def clean_adjudicated_category(row):
            if row["Primary Adjudicated Category"] != "NA":
                return row["Primary Adjudicated Category"]
            else:
                # The row was empty, so perhaps this sample has not yet been adjudicated.
                # However, we also need to check that this was not a "Healthy Reference" sample,
                # as these never go through the adjudication process.
                if row["EnrollmentCategory"] in ["Healthy Reference"]:
                    return "Healthy Reference"
                return ""
        
        # Note: we clean enrollment category first, as the adjudicated category cleaning depends on the cleaned enrollment category values.
        adata.obs["EnrollmentCategory"] = adata.obs.apply(clean_enrollment_category, axis='columns')
        adata.obs["AdjudicatedCategory"] = adata.obs.apply(clean_adjudicated_category, axis='columns')
        

        # TODO: expand all acronyms in the sample group names, to be consistent?

        # When using the adjudicated sample categorizations,
        # we want to aggregate the AKI subcategories of ATI and AIN (roll them up into a parent "AKI" category),
        # since there are currently not enough samples for comparisons at the subcategory level to make sense.
        # To do so, we add a new "MergedAdjudicatedCategory" column that merges ATI and AIN into AKI.
        def merge_ati_and_ain(row):
            if row["AdjudicatedCategory"] in ["Acute Tubular Injury", "Acute Interstitial Nephritis"]:
                return "AKI"
            else:
                return row["AdjudicatedCategory"]
        adata.obs["MergedAdjudicatedCategory"] = adata.obs.apply(merge_ati_and_ain, axis='columns')

        # TODO: process other clinical columns? Sex, age group, etc.

        adata.obs = adata.obs.rename(columns={
            "subclass.l1": "subclass_l1",
            "subclass.l2": "subclass_l2",
            "SpecimenID": "specimen",
            "KPMP Participant ID": "patient",
        })

        for colname in adata.obs.columns:
            if pd.api.types.is_string_dtype(adata.obs[colname]) or str(adata.obs[colname].dtype) == "object":
                print(f"Filling NAs in string column {colname} with 'NA'")
                adata.obs[colname] = adata.obs[colname].fillna("NA")
            else:
                print(f"Not filling NAs in non-string column {colname} of type {adata.obs[colname].dtype}")

        # Column names cannot contain slashes
        adata.obs = adata.obs.rename(columns=dict(zip(adata.obs.columns, [c.replace("/", " per ") for c in adata.obs.columns])))





        # TODO: Use dask to convert the scipy.sparse array to a dense numpy array, to avoid this memory spike?
        # TODO: Also use LazyAnnData.put_da_to_zarr_layer to save the dense array?
        # # Perhaps save "counts" to "counts_sparse" layer first, then convert to dense, then save result to "counts" layer.
        # dask_sparse_array = da.from_array(large_sparse_array, chunks=(1000, 1000))
        # # Define a function to convert a sparse chunk to a dense NumPy array
        # def sparse_to_dense_chunk(chunk):
        #     return chunk.toarray()
        # # Apply the function to each block using map_blocks
        # # The 'meta' argument helps Dask infer the output type and shape
        # dask_dense_array = dask_sparse_array.map_blocks(
        #     sparse_to_dense_chunk,
        #     dtype=large_sparse_array.dtype,
        #     meta=np.array([]) # Dask expects a NumPy array for meta
        # )

        # Reference: https://stackoverflow.com/questions/30416695/numpy-and-scipy-difference-between-todense-and-toarray

        # DO NOT CONVERT TO DENSE, KEEP AS SPARSE.
        # adata.layers["counts_dense"] = adata.layers["counts"].toarray()

        return adata

    # TODO: get this info from the Snakemake config yaml file
    donor_id_col = "patient"
    sample_id_col = "specimen"
    sample_group_pairs = [
        # AKI vs. HRT
        ('EnrollmentCategory', ('Healthy Reference', 'AKI')),
        # AKI vs. H-CKD. (H-CKD not in enrollment category values anymore. Should I use "Hypertension History" Yes/No column?)
        ('EnrollmentCategory', ('AKI', 'CKD')),
        # D-CKD vs. HRT. (D-CKD not in enrollment category values anymore. Should I use "Diabetes History" Yes/No column?)
        ('EnrollmentCategory', ('CKD', 'Healthy Reference')),
        # DM-R comparisons
        ('EnrollmentCategory', ('Healthy Reference', 'DM-R')),
        ('EnrollmentCategory', ('CKD', 'DM-R')),
        ('EnrollmentCategory', ('AKI', 'DM-R')),
        # Diabetes CKD vs. Hypertension CKD. (DKD nor H-CKD not in enrollment category values anymore. Should I use Yes/No columns?)
        #('EnrollmentCategory', ('DKD', 'H-CKD')),
        # D-CKD vs. HRT
        ('AdjudicatedCategory', ('Diabetic Kidney Disease', 'Healthy Reference')),
        # Acute tubular injury vs. HRT
        ('AdjudicatedCategory', ('Acute Tubular Injury', 'Healthy Reference')),
        # Acute interstitial nephritis vs. HRT
        ('AdjudicatedCategory', ('Acute Interstitial Nephritis', 'Healthy Reference')),
        # Diabetes CKD vs. Hypertension CKD
        ('AdjudicatedCategory', ('Diabetic Kidney Disease', 'Hypertensive Kidney Disease')),
        # ATN vs. AIN
        ('AdjudicatedCategory', ('Acute Interstitial Nephritis', 'Acute Tubular Injury')),

        # Merged AKI (ATI and AIN) vs. others
        ('MergedAdjudicatedCategory', ('AKI', 'Healthy Reference')),
        ('MergedAdjudicatedCategory', ('AKI', 'Diabetic Kidney Disease')),
        ('MergedAdjudicatedCategory', ('AKI', 'Hypertensive Kidney Disease')),

        # TODO: use Diabetes History and Hypertension History columns here.

        # TODO: For healthy: Define as HRT samples from KPMP AND (either kidney stone or pre-perfusion biopsies).
        # I.e., exclude tumor nephrectomies from healthy. Exclude non-KPMP biopsies.
    ]
    cell_type_cols = [
        "subclass_l2",
        "subclass_l1",
    ]

    run_all(
        get_adata,
        out_h5ad_path=args.output,
        out_zarr_path=args.output_zarr,
        overwrite=args.overwrite,
        client=create_o2_dask_client(memory_limit=args.mem_limit),
        donor_id_col=donor_id_col,
        sample_id_col=sample_id_col,
        sample_group_pairs=sample_group_pairs,
        cell_type_cols=cell_type_cols,
        stop_early=args.stop_early,
    )

    print("Done")
