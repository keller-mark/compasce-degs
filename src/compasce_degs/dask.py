import os


def create_dask_client(dask_temp_dir=None, n_workers=2, threads_per_worker=2, memory_limit='24GB'):
    import dask
    from dask.distributed import Client, LocalCluster

    if dask_temp_dir is not None:
        os.makedirs(dask_temp_dir, exist_ok=True)
        dask.config.set({ "temporary_directory": dask_temp_dir })
    
    dask.config.set({ 'logging.distributed.scheduler': 'error' })

    # Should request at least 96GB of memory for this job.
    # TODO: more kwargs for LocalCluster? option to pass a cluster instead?
    cluster = LocalCluster(n_workers=n_workers, threads_per_worker=threads_per_worker, memory_limit=memory_limit)
    client = Client(cluster)

    return client


def create_dask_wrapper(inner_func):
    def dask_wrapper(get_input_arr, put_output_arr, **kwargs):
        arr = get_input_arr() # returns dask arr
        output_arr = inner_func(arr, **kwargs)
        return put_output_arr(output_arr) # accepts dask arr
    return dask_wrapper
