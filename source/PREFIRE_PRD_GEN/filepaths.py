# This is a centralized module (within this package) where other modules can
#  obtain selected path information and default filepaths.

import os.path

code_dir = os.path.dirname(os.path.realpath(__file__))
package_dir = os.path.abspath(os.path.join(code_dir, "..", ".."))
package_ancillary_data_dir = os.path.abspath(os.path.join(package_dir,
                                                           "dist", "ancillary"))

scipkg_version_fpath = os.path.join(package_dir, "VERSION.txt")
scipkg_prdgitv_fpath = os.path.join(package_dir, "dist", "prdgit_version.txt")

_shared_geometry_filespec_fpath = os.path.join(package_ancillary_data_dir,
                                              "shared_geometry_filespecs.json")
_final_products_filespec_fpath = os.path.join(package_ancillary_data_dir,
                                               "final_products_filespecs.json")
_product_ops_to_perform_fpath = os.path.join(package_ancillary_data_dir,
                                                 "product_ops_to_perform.json")
_extracted_metrics_filespec_fpath = os.path.join(package_ancillary_data_dir,
                                            "extracted_metrics_filespecs.json")
