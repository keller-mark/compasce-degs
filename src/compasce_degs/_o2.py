import os
import platform

from .dask import create_dask_client


def is_cluster():
    is_gh_actions = os.environ.get("CI", False)
    is_macbook = (platform.system() == "Darwin")
    return not is_gh_actions and not is_macbook


def create_o2_dask_client(**kwargs):
    dask_temp_dir = None
    if is_cluster():
        O2_USER = os.environ["USER"]
        dask_temp_dir = f"/n/scratch/users/{O2_USER[0]}/{O2_USER}/scmd-analysis-tmp"
    return create_dask_client(dask_temp_dir=dask_temp_dir, **kwargs)
