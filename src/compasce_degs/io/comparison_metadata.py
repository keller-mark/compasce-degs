import json
import numpy as np
import zarr
from os.path import join
from .cdata import dir_name_to_str


# Reference: https://stackoverflow.com/a/57915246
class NpEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return super(NpEncoder, self).default(obj)

# Reference: https://stackoverflow.com/a/68557484
def deep_update(mapping, *updating_mappings) :
    updated_mapping = mapping.copy()
    for updating_mapping in updating_mappings:
        for k, v in updating_mapping.items():
            if k in updated_mapping and isinstance(updated_mapping[k], dict) and isinstance(v, dict):
                updated_mapping[k] = deep_update(updated_mapping[k], v)
            elif k in updated_mapping and isinstance(updated_mapping[k], list) and isinstance(v, list):
                updated_mapping[k] = [
                    *updated_mapping[k],
                    *v
                ]
            else:
                updated_mapping[k] = v
    return updated_mapping


# See https://observablehq.com/d/e7c03bf319f20f86
class ComparisonMetadata:
    def __init__(self, comparison_key):
        self.comparison_key = comparison_key
        self.comparison_key_str = dir_name_to_str(comparison_key)
        self.items = []
    
    def get_df_key(self, df_type):
        return f"{self.comparison_key_str}.{df_type}"
    
    def append_df(self, adata_key, df_type, df_params, df_c_vals):
        self.items.append({
            "path": join(adata_key, self.get_df_key(df_type)),
            "coordination_values": df_c_vals,
            "analysis_type": df_type,
            "analysis_params": df_params,
        })
        return self.get_df_key(df_type)
    
    def get_dict(self):
        return {
            self.comparison_key_str: {
                "comparison": self.comparison_key,
                "results": self.items,
            }
        }
    
class MultiComparisonMetadata:
    def __init__(self, sample_group_pairs=None, sample_id_col=None, donor_id_col=None, cell_type_cols=None):
        self.schema_version = "0.0.2"
        self._comparisons = []
        self.sample_group_pairs = sample_group_pairs
        self.sample_id_col = sample_id_col
        self.donor_id_col = donor_id_col
        self.cell_type_cols = cell_type_cols
        self._prev_comparisons_dict = dict()
    
    def load_state(self, zarr_path, include_comparisons=True):
        z = zarr.open(zarr_path, mode="a")
        if "/uns/comparison_metadata" in z:
            prev = json.loads(str(z["/uns/comparison_metadata"][()]))
            if include_comparisons:
                self._prev_comparisons_dict = prev["comparisons"]
            if self.sample_group_pairs is None:
                self.sample_group_pairs = prev["sample_group_pairs"]
            if self.sample_id_col is None:
                self.sample_id_col = prev["sample_id_col"]
            if self.donor_id_col is None:
                self.donor_id_col = prev["donor_id_col"]
            if self.cell_type_cols is None:
                self.cell_type_cols = prev["cell_type_cols"]
    
    def merge_states(self, zarr_path, suffixes=None):
        z = zarr.open(zarr_path, mode="a")
        comparisons_dict = self._prev_comparisons_dict
        for suffix in suffixes:
            if f"/uns/comparison_metadata.{suffix}" in z:
                prev = json.loads(str(z[f"/uns/comparison_metadata.{suffix}"][()]))
                comparisons_dict = deep_update(comparisons_dict, prev["comparisons"])
        self._prev_comparisons_dict = comparisons_dict
    
    def add_comparison(self, comparison_key):
        c = ComparisonMetadata(comparison_key)
        self._comparisons.append(c)
        return c
    
    def serialize(self):
        comparisons_dict = self._prev_comparisons_dict
        for c in self._comparisons:
            comparisons_dict = deep_update(comparisons_dict, c.get_dict())
        return json.dumps({
            "schema_version": self.schema_version,
            "comparisons": comparisons_dict,
            "sample_id_col": self.sample_id_col,
            "sample_group_pairs": self.sample_group_pairs,
            "donor_id_col": self.donor_id_col,
            "cell_type_cols": self.cell_type_cols,
        })
    