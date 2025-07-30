#!/usr/bin/env bash

## IMPORTANT: Only run this script from the directory it resides in, i.e. with
##             ./run.sh    OR    bash run.sh

##===========================================================================##
## This script contains hardwired information necessary for this algorithm's
##  delivery to and testing within the SDPS (Science Data Processing System).
##
## ** In general, do not push changes to this file to its primary git
##     repository (exceptions include adding a new environment var for
##     algorithm config) **
##
## ++ Instead, make a LOCAL copy of this script (e.g., my_run.sh; do not
##     push that local copy to the primary git repository either) and modify
##     and run that for general algorithm testing and development.
##===========================================================================##

absfpath() {
  # Generate absolute filepath from a relative (or even an absolute) filepath.
  #
  # Based on (circa Oct 2023) https://stackoverflow.com/questions/3915040/how-to-obtain-the-absolute-path-of-a-file-via-shell-bash-zsh-sh
  # 
  # $1     : a relative (or even an absolute) filepath
  # Returns the corresponding absolute filepath.
  if [ -d "$1" ]; then
    # dir
    (cd "$1"; pwd)
  elif [ -f "$1" ]; then
    # file
    if [[ $1 = /* ]]; then
      echo "$1"
    elif [[ $1 == */* ]]; then
      echo "$(cd "${1%/*}"; pwd)/${1##*/}"
    else
      echo "$(pwd)/$1"
    fi
  fi
}

activate_conda_env () {
  . "$1"/bin/activate;
}

deactivate_conda_env () {
  . "$1"/bin/deactivate;
}

#set -ve;  # Exit on the first error, and print out commands as we execute them
set -e;  # Exit on the first error

# Determine the absolute path of the current working directory:
#  (this is typically the package test/ directory)
readonly base_dir="$(absfpath ".")";

hn=`hostname -s`;  # Hostname

# NOTE: Set the input/output directories to absolute paths (relative to the
#        current working directory, 'base_dir').

non_SDPS_hostname="longwave";

srcf_dir="${base_dir}/inputs";

cfg_str1="01-${srcf_dir}/prefire_01_payload_tlm_2024_07_25_18_33_27.bin";
cfg_str2="${srcf_dir}/raw-PREFIRE_SAT2_1B-RAD_P01_R00_20241009180240_02077.nc";
cfg_str3="${srcf_dir}/raw-PREFIRE_SAT2_2B-MSK_P00_R00_20241007071543_02040.nc";

# Specify that numpy, scipy, et cetera should not use more than one thread or
#  process):
MKL_NUM_THREADS=1;
NUMEXPR_NUM_THREADS=1;
OMP_NUM_THREADS=1;
VECLIB_MAXIMUM_THREADS=1;
OPENBLAS_NUM_THREADS=1;
export MKL_NUM_THREADS NUMEXPR_NUM_THREADS OMP_NUM_THREADS;
export VECLIB_MAXIMUM_THREADS OPENBLAS_NUM_THREADS;

# Some environment vars that convey configuration info to the algorithm:

this_top_dir="$(absfpath "${base_dir}/..")";

PACKAGE_TOP_DIR="${this_top_dir}";
ANCILLARY_DATA_DIR="${this_top_dir}/dist/ancillary";

OUTPUT_PRODUCT_DIR="${base_dir}/outputs";

export PACKAGE_TOP_DIR ANCILLARY_DATA_DIR OUTPUT_PRODUCT_DIR;

  # syntax is:   key1:value(s) | key2:value(s) | ...
  #   Examples of key:value combinations:
  #      "lat_bounds_for_asc_pass: -90. to -60., 60. to 90."
  #      "lat_bounds_for_desc_pass: -90. to -60., 60. to 90."
  #      "lat_bounds_for_all_passes: -90. to -60., 60. to 90."
OUTPUT_OPTIONS="";
export OUTPUT_OPTIONS;

# Check if output file directory exists; if not, bail:
tmpdir=${OUTPUT_PRODUCT_DIR};
test -d "${tmpdir}" || { echo "Output directory does not exist: ${tmpdir}"; exit 1; }

# If custom conda environment files exist, activate that conda environment:
conda_env_dir="${this_top_dir}/dist/c_env_for_PREFIRE_PRD_GEN";
if [ -d "${conda_env_dir}" ]; then
   activate_conda_env "${conda_env_dir}";
fi

# Execute script that writes a new 'prdgit_version.txt', which contains
#  product moniker(s) and current (latest) git hash(es) that are part of the
#  provenance of this package's product(s).
# *** This step should not be done within the SDPS, since that file is
#     created just before delivery to the SDPS.
if [ ! -f "${this_top_dir}/dist/for_SDPS_delivery.txt" ]; then
   python "${this_top_dir}/dist/determine_prdgit.py";
fi

for cfg_str in ${cfg_str1} ${cfg_str2} ${cfg_str3}
do
   # * Will concatenate into a single file if more than one file matches
   #    INPUT_FPATH_OR_FPREFIX.
   # * If there is only a single match, or after concatenation, then produce a
   #    final product file (and any attendant metadata files) -- but only if
   #    the product type is found in this package's finalized-product filespecs.
   INPUT_FPATH_OR_FPREFIX="${cfg_str}";

   export INPUT_FPATH_OR_FPREFIX;

   # Execute primary driver:
   python "${this_top_dir}/dist/produce_PRD.py";
done

# If custom conda environment files exist, DEactivate that conda environment:
if [ -d "${conda_env_dir}" ]; then
   deactivate_conda_env "${conda_env_dir}";
fi

