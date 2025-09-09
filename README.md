# PREFIRE_PRD_GEN (Product Generator)

Python-based package for vetting/preparing final NetCDF-format PREFIRE AUX, Level-1A, Level-1B, Level-2B, and Level-3 output products, generating their (and Level-0) attendant archival metadata files, representative-visualization images for NASA Earthdata search, selected-visualization COG (cloud-optimized GeoTIFF) images for NASA Worldview, and extracting metrics for use in trending the hardware/radiometric/algorithm status. It also can simply concatenate multiple along-track file chunks (NetCDF-format) into a single file that contains a full granule's worth of those file contents (even for non-final-product intermediate files).

This package is intended to be run immediately after each full granule's PREFIRE science package processing job ensemble (i.e., multiple concurrent and/or sequential jobs, each processing part of the same granule) is finished.  This package may simply construct an intermediate product file (from one or more separate parts) or construct a final product file (from one or more separate parts).

Many of the modules in this package are also intended to be used by any Python code (where appropriate) that may be a part of each of the science product algorithm packages, in order to reduce code duplication when reading, writing, et cetera.

This code is released under the terms of this [LICENSE](LICENSE).  The version of this package can be found in [VERSION.txt](VERSION.txt).

# Installation

## Requirements

Python version 3.8+ is required, along with the following third-party Python packages: numpy, netcdf4, matplotlib, cartopy, shapely, geopandas, geocube, skyfield

The associated (Python-based) git repository ['PREFIRE_tools'](https://github.com/UW-PREFIRE/PREFIRE_tools) is also required for the proper operation of this package.

## Python Environment Setup

It is recommended to install the above Python packages in a dedicated conda environment (or something similar).  The packages used (and their versions) can be found in [conda_env.list](dist/conda_env.list).

For example, using conda (and specifying Python 3.10.x from the conda-forge channel):

```
conda create --name for_PREFIRE_PRD_GEN -c conda-forge python=3.10;
conda activate for_PREFIRE_PRD_GEN;
conda install -c conda-forge numpy netcdf4 matplotlib cartopy shapely geopandas geocube;
```

The location of 'PREFIRE_tools' depends on the value of the user's PYTHONPATH and/or sys.path -- for example, one could simply add each of those git repositories' local root Python source code directory to PYTHONPATH.

Operationally, however, this package uses symbolic links to those git repositories' local root Python source code directories (or full copies of the same) in the source/ directory.  To use this symlink method (assuming that all PREFIRE code repositories are in the same parent directory, and that the PYTHONPATH environment variable is unset or empty):

```
cd source;
ln -s ../../PREFIRE_tools/source/python/PREFIRE_tools PREFIRE_tools;
```

## Environment Variables

### Each job (executing this science algorithm package) is configured via information contained within environment variables.

### To specify that numpy, geopandas, et cetera used by this algorithm should not use more than one thread or process, the below environment variables are expected to be set:

```
MKL_NUM_THREADS=1
NUMEXPR_NUM_THREADS=1
OMP_NUM_THREADS=1
VECLIB_MAXIMUM_THREADS=1
OPENBLAS_NUM_THREADS=1
```

### The following environment variables configure the run (also see test/run.sh or test/run.ps1):

ANCILLARY_DATA_DIR  :  the package's ancillary data directory (should be an absolute path)
   
INPUT_FPATH_OR_PREFIX  :  the absolute filepath of the input granule; _OR_ the filepath prefix (absolute path plus the first part of the filename) of the input granule part files; _OR_ (for non-NetCDF-format files only) "xx-absolute_filepath", where 'absolute_filepath' is that of the input granule, and 'xx' is a numeric revision/collection version (starts at 01, next is 02, and so on...)

OUTPUT_PRODUCT_DIR  :  the directory in which all meaningful output will be written (should be an absolute path)

OUTPUT_OPTIONS  :  (optional) options that control the output.  The general syntax is "key1:value(s) | key2:value(s) | ..." -- Examples of key:value combinations: (1) "lat_bounds_for_asc_pass: -90. to -60., 60. to 90.", (2) "lat_bounds_for_desc_pass: -90. to -60., 60. to 90.", (3) "lat_bounds_for_all_passes: -90. to -60., 60. to 90."

# Running the test script(s)

## Obtain and unpack any ancillary data

None (for this version).

## Prepare the output directory:

`cd test;`

On Linux/UNIX systems, possibly create a useful symbolic link to the test input data (if needed):

`ln -s WHEREEVER_THE_DATA_IS/inputs inputs;`

Prepare the output directory (Linux/UNIX example):

`mkdir -p outputs;`

_OR_ perhaps something like

`ln -s /data/users/myuser/data-PREFIRE_PRD_GEN/outputs outputs;`

## Run the PRD_GEN package

### A Linux/UNIX example

`cp run.sh my-run.sh;`

Edit `my-run.sh` as needed (e.g., change input file names)

`./my-run.sh`

The output file(s) will be in `test/outputs/`

### _The creation of this code was supported by NASA, as part of the PREFIRE (Polar Radiant Energy in the Far-InfraRed Experiment) CubeSat mission._
