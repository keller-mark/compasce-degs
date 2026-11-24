from os.path import join, dirname
import platform

def join_arr(*args):
  return join(*args, ".zarray")

def join_attrs(*args):
  return join(*args, ".zattrs")

def join_group(*args):
  return join(*args, ".zgroup")

def join_zdone(*args):
  return join(*args, ".zdone")

def resolve_zarr(p):
    assert ".zarr" in p
    if p.endswith(".zarr"):
        return p
    else:
        return resolve_zarr(dirname(p))


# Check if this is running on O2
IS_O2 = (platform.system() == "Linux")

# Directory / file constants
SRC_DIR = "src"
SCRIPTS_DIR = "scripts"
DATA_DIR = ("data" if not IS_O2 else "/n/data1/hms/dbmi/gehlenborg/lab/scmd-analysis")
RAW_DIR = join(DATA_DIR, "raw")
COMMON_INTERMEDIATE_DIR = join(DATA_DIR, "intermediate")
COMMON_PROCESSED_DIR = join(DATA_DIR, "processed")