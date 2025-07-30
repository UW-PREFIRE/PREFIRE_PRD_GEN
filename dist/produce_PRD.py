"""
Produce final NetCDF-format PREFIRE AUX, Level-1A, Level-1B, Level-2B, and
 Level-3 output products, as well as generating their attendant archival
 metadata files.  Also can simply concatenate multiple along-track file chunks
 into a single file that contains a full granule's worth of those file contents
 (even for non-final-product intermediate files).

This program requires Python version 3.6 or later, and is importable as a 
python module.
"""

  # From the Python standard library:
from pathlib import Path
import os
import sys
import argparse

  # From other external Python packages:

  # Custom utilities:


#--------------------------------------------------------------------------
def main():
    """Driver routine."""

    package_top_Path = Path(os.environ["PACKAGE_TOP_DIR"])

    sys.path.append(str(package_top_Path / "source"))
    from PREFIRE_PRD_GEN.apply_PRD_GEN import apply_PRD_GEN

    input_fpath_or_fpathpfx = os.environ["INPUT_FPATH_OR_FPREFIX"]
    output_Path = Path(os.environ["OUTPUT_PRODUCT_DIR"])
    ancillary_Path = Path(os.environ["ANCILLARY_DATA_DIR"])

    output_options = None
    if "OUTPUT_OPTIONS" in os.environ:
        tmp = os.environ["OUTPUT_OPTIONS"].strip()
        if len(tmp) > 0:
            output_options = tmp

    apply_PRD_GEN(input_fpath_or_fpathpfx, output_Path, ancillary_Path,
                  package_top_Path, output_options)


if __name__ == "__main__":
    # Process arguments:
    arg_description = ("Produce concatenated and/or final forms of PREFIRE "
                       "science package output files.")
    arg_parser = argparse.ArgumentParser(description=arg_description)

    args = arg_parser.parse_args()

    # Run driver:
    main()
