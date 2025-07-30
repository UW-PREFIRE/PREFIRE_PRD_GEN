## IMPORTANT: Only run this script from the directory it resides in, i.e. with
##             .\run.ps1

##===========================================================================##
## This script contains hardwired information necessary for this algorithm's
##  delivery to and testing within the SDPS (Science Data Processing System).
##
## ** In general, do not push changes to this file to its primary git
##     repository (exceptions include adding a new environment var for
##     algorithm config) **
##
## ++ Instead, make a LOCAL copy of this script (e.g., my_run.ps1; do not
##     push that local copy to the primary git repository either) and modify
##     and run that for general algorithm testing and development.
##===========================================================================##

# Determine the absolute path of the current working directory:
#  (this is typically the package test/ directory)
$base_dir = $pwd

# NOTE: Set the input/output directories to absolute paths (relative to the
#        current working directory, 'base_dir').

$srcf_dir = "$base_dir\inputs"

$cfg_str1 = "01-$srcf_dir\prefire_01_payload_tlm_2024_07_25_18_33_27.bin"
$cfg_str2 = "$srcf_dir\raw-PREFIRE_SAT2_1B-RAD_P01_R00_20241009180240_02077.nc"
$cfg_str3 = "$srcf_dir\PREFIRE_SAT2_2B-MSK_P00_R00_20241007040521_02038.nc"

# Specify that numpy, scipy, et cetera should not use more than one thread or
#  process):
$env:MKL_NUM_THREADS = '1'
$env:NUMEXPR_NUM_THREADS = '1'
$env:OMP_NUM_THREADS = '1'
$env:VECLIB_MAXIMUM_THREADS = '1'
$env:OPENBLAS_NUM_THREADS = '1'

# Some environment vars that convey configuration info to the algorithm:

$this_top_dir = [IO.Path]::GetFullPath("$base_dir\..")

$env:PACKAGE_TOP_DIR = "$this_top_dir"
$env:ANCILLARY_DATA_DIR = "$this_top_dir\dist\ancillary"

$env:OUTPUT_PRODUCT_DIR = "$base_dir\outputs"

  # syntax is:   key1:value(s) | key2:value(s) | ...
  #   Examples of key:value combinations:
  #      "lat_bounds_for_asc_pass: -90. to -60., 60. to 90."
  #      "lat_bounds_for_desc_pass: -90. to -60., 60. to 90."
  #      "lat_bounds_for_all_passes: -90. to -60., 60. to 90."
#$env:OUTPUT_OPTIONS = ""

# Check if output file directory exists; if not, bail:
$tmpdir = "$env:OUTPUT_PRODUCT_DIR"
If (-not (Test-Path -Path $tmpdir)) {
  throw "Output directory does not exist: $tmpdir"
}

# Execute script that writes a new 'prdgit_version.txt', which contains
#  product moniker(s) and current (latest) git hash(es) that are part of the
#  provenance of this package's product(s).
# *** This step should not be done within the SDPS, since that file is
#     created just before delivery to the SDPS.
If (-not (Test-Path -Path "$this_top_dir\dist\for_SDPS_delivery.txt")) {
  python "$this_top_dir\dist\determine_prdgit.py"
}

$items = @($cfg_str1, $cfg_str2, $cfg_str3)
foreach ($cfg_str in $items) {
  # * Will concatenate into a single file if more than one file matches
  #    INPUT_FPATH_OR_FPREFIX.
  # * If there is only a single match, or after concatenation, then produce a
  #    final product file (and any attendant metadata files) -- but only if
  #    the product type is found in this package's finalized-product filespecs
  $env:INPUT_FPATH_OR_FPREFIX = $cfg_str

  # Execute primary driver:
  python "$this_top_dir\dist\produce_PRD.py"
}
