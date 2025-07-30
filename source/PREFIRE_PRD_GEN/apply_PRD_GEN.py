"""
This program requires python version 3.6 or later, and is importable as a
python module.
"""

  # From the Python standard library:
from pathlib import Path
import glob
import json
import datetime
import os
import re
from itertools import permutations
import copy

  # From other external Python packages:
import netCDF4

  # Custom utilities:
from PREFIRE_PRD_GEN.file_concatenation import concatenate_along_dim
from PREFIRE_PRD_GEN.file_description import produce_DAAC_metadata
from PREFIRE_PRD_GEN.file_finalization import finalize_product_file
from PREFIRE_PRD_GEN.metrics_extraction import extract_metrics
import PREFIRE_PRD_GEN.filepaths as filepaths


non_nc_info = ({"nospecs_TIRS_L0_Data":
                  r"prefire_0\d_payload_tlm_[0-9]+_[0-9]+_[0-9]+_[0-9]+_[0-9]+_[0-9]+\.bin",
                "nospecs_TIRS_L0_CrtdData":
                  r"prefire_0\d_payload_tlm_[0-9]+_[0-9]+_[0-9]+\.bin",
                "nospecs_SCbus_L0_CrtdData":
                  r"prefire_0\d_bus_tlm_[0-9]+_[0-9]+_[0-9]+\.csv"})


#--------------------------------------------------------------------------
def input_file_processing_plan(input_fpath_l):
    """
    Determine and return a processing plan for the input file(s), taking into
    account file type, content type, etc.
    """

    # Initial/default values:
    input_is_all_nc = False
    in_file_type = {}
    processing_plan = {}
    in_file_info = {}

    # Read in product-operations-to_perform specification info:
    specs_json_fpath = filepaths._product_ops_to_perform_fpath
    with open(specs_json_fpath, 'r') as f:
        tmp = json.load(f)
    product_ops_to_perform = copy.deepcopy(tmp)
    for key in tmp:
        tokens = tuple([x.strip() for x in key.split(',')])
        if len(tokens) > 1:
            for perm in permutations(tokens[1:]):
                for i in range(len(tokens)):
                    p_key_t = perm[:i] + tokens[0:1] + perm[i:]
                    p_key = ", ".join(p_key_t)
                    product_ops_to_perform[p_key] = product_ops_to_perform[key]

    # All NetCDF-format input?
    try:
        for in_fpath in input_fpath_l:
            with netCDF4.Dataset(in_fpath, 'r') as nc_ds:
                dims = nc_ds.dimensions
        input_is_all_nc = True
        in_file_type["NetCDF-format file"] = True
    except:
        in_file_type["other-format file"] = True

    # Is the amount of input files allowable?
    if len(input_fpath_l) > 1:
        if input_is_all_nc:
            processing_plan["concatenate"] = True
        else:
            msg = ("One or more other-format input files found, but this "
                   "program can only concatenate NetCDF-format files (list "
                   "of all input files: {})".format(", ".join(input_fpath_l)))
            raise RuntimeError(msg)

    # Get product type key (setting other info along the way) for NetCDF-format
    #  input file(s):
    if input_is_all_nc:
          # Read in final-product file specification info:
        filespecs_json_fpath = filepaths._final_products_filespec_fpath
        in_file_info["filespecs_json_fpath"] = filespecs_json_fpath
        with open(filespecs_json_fpath, 'r') as f:
            filespecs_of_final_products = json.load(f)

        # Are the non-geometry groups in the input file(s) all in the filespecs
        #  for finalized files?
        unfinalized_type = None  # Default
        product_type_l = []
        ig_save_l = []
        in_file_info["use_shared_geom"] = False
        with netCDF4.Dataset(input_fpath_l[0], 'r') as nc_ds:
            in_file_info["SV_IDstr"] = nc_ds.spacecraft_ID[-1]
            for inpf_gname in nc_ds.groups:
                if ("Attributes" not in inpf_gname) and (
                       "Geometry" not in inpf_gname):  # i.e., a "normal" group
                    ig_save_l.append(inpf_gname)
                    for outf_gname in filespecs_of_final_products:
                        if inpf_gname == outf_gname:
                            product_type_l.append(outf_gname)
                    if (inpf_gname == "TIRS_L0_Data" or
                        inpf_gname == "SCbus_L0_Data" or
                        inpf_gname == "Orbit_L0_Data"):
                            unfinalized_type = inpf_gname
                else:  # Attribute group or Geometry group
                    if (inpf_gname == "Geometry" or
                                                inpf_gname == "ProtoGeometry"):
                        in_file_info["use_shared_geom"] = True
        in_file_info["ig_save_l"] = ig_save_l

        if len(product_type_l) > 0:
            processing_plan["finalize"] = True
            in_file_info["product_type_l"] = product_type_l
            product_type_key = ", ".join(product_type_l)
        elif unfinalized_type is not None:
            product_type_key = unfinalized_type
        else:
            msg = ("Cannot finalize product file because the "
                   "existing groups ({}) were not found in the filespecs "
                  "file (dist/ancillary/*.json).".format(", ".join(ig_save_l)))
            raise RuntimeError(msg)

    # Get product type key (setting other info along the way) for an
    #  other-format input file:
    else:
        input_fn = os.path.basename(input_fpath_l[0])
        found = False
        for key in non_nc_info:
            id_pattern = non_nc_info[key]
            pattern = re.compile(id_pattern)
            if pattern.match(input_fn):
                found = True
                in_file_info["SV_IDstr"] = input_fn[9]
                if "CrtdData" in key:
                    tmp_l = [input_fn[-48:-35+1], input_fn[-33:-20+1]]
                    for k, s in zip(["UTC_c_start", "UTC_c_end"], tmp_l):
                        in_file_info[k] = ("{}-{}-{}T{}:{}:{}.000000".format(
                                                  s[0:4], s[4:6], s[6:8],
                                                  s[8:10], s[10:12], s[12:14]))
                elif key == "nospecs_TIRS_L0_Data":
                    # Crack open all CCSDS packets and find this product file's
                    #  UTC start and end:
                    from PREFIRE_tools.utils.time import ctime_to_UTC_DT, \
                                                  init_leap_s_for_ctimeRefEpoch
                    import PREFIRE_tools.TIRS.TIRS_packet as TIRSp
                    import PREFIRE_tools.utils.CCSDS_packet_header as Cph

                    leap_s_info = init_leap_s_for_ctimeRefEpoch(
                                            [2000, 1, 1, 0, 0 ,0],
                                            epoch_for_ctime_is_actual_UTC=True)
                    bytesize_of_CCSDS_hdr = Cph.CCSDS_header_specs()
                    payload_scipkt_buffinfo = (
                                         TIRSp.get_scipkt_bytebuffer_indices(
                                                        bytesize_of_CCSDS_hdr))
                    nonscipkt_log = []
                    ctime_s_l = []

                    in_f = open(input_fpath_l[0], "rb")
                    scipkt_cnt = -1
                    while True:
                        scipkt_data_t = TIRSp.read_one_scipkt(in_f,
                                bytesize_of_CCSDS_hdr, payload_scipkt_buffinfo,
                                                    nonscipkt_log, leap_s_info)
                        if scipkt_data_t[0] is None:
                            break

                        _, _, encoder_pos, ctime_s, _ = scipkt_data_t
                        if ctime_s > 6.97e8:  # After approximately 2022-02-01
                            ctime_s_l.append(ctime_s)

                        scipkt_cnt += 1
                    in_f.close()

                    UTC_DT, _ = ctime_to_UTC_DT(
                                              [min(ctime_s_l), max(ctime_s_l)],
                                                's', leap_s_info)
                    in_file_info["UTC_c_start"] = (
                                    UTC_DT[0].strftime("%Y-%m-%dT%H:%M:%S.%f"))
                    in_file_info["UTC_c_end"] = (
                                    UTC_DT[1].strftime("%Y-%m-%dT%H:%M:%S.%f"))
                else:
                    in_file_info["UTC_c_start"] = None  # No UTC coverage bounds
                    in_file_info["UTC_c_end"] = None    #  available
                product_type_key = key
                break
        if not found:
            raise RuntimeError(
                             f"Unknown type of input file {input_fpath_l[0]}")

    in_file_type["product_type_key"] = product_type_key

    # Determine product moniker and the operations to perform for it:
    try:
        product_moniker = (
                   product_ops_to_perform[product_type_key]["product_moniker"])
        product_operations = (
                        product_ops_to_perform[product_type_key]["operations"])
    except KeyError:
        msg = ("ERROR: The product type is unknown for the existing input "
               "product key(s) ({})".format(product_type_key))
        raise RuntimeError(msg)
      # Special cases:
    if product_moniker == "1B-GEOM":
        if "1B-NTCG" in os.path.basename(input_fpath_l[0]):
            product_moniker = "1B-NTCG"
    elif product_moniker == "1B-RAD":
        if "1B-NLRAD" in os.path.basename(input_fpath_l[0]):
                product_moniker = "1B-NLRAD"
                create_DAAC_mdf = False
    in_file_type["product_moniker"] = product_moniker
    for op in product_operations:
        processing_plan[op] = True

    return (in_file_type, processing_plan, in_file_info)


#--------------------------------------------------------------------------
def apply_PRD_GEN(input_fpath_or_fpathpfx, output_Path, ancillary_Path,
                  package_top_Path, output_options):
    """
    Applies some Product Generator operations to input file(s) -- NetCDF-format
     concatenation, finalization, extraction of trending-metrics, and/or
     production of DAAC per-granule metadata, representative-visualization
     images for NASA Earthdata search, selected-visualization COG
     (cloud-optimized GeoTIFF) images for NASA Worldview -- and produces one or
     more new files as output.

    Input Parameters
    ----------
    input_fpath_or_fpathpfx : str
        Input absolute filepath (any file format) or the common absolute
         filepath prefix of the NetCDF-format chunkfiles to be concatenated.
        Note: The archival version ID (two-digit; e.g. "01") can be specified
         for non-NetCDF-format files by putting it at the beginning of this
         string, delimited by a '-' (e.g., "01-prefire_02_payload_...").  In
         this case, those first 3 characters will be stripped from the
         fpath/fpath_pfx before use.
    output_Path : pathlib.Path object
        Directory in which to put the output file(s). 
    ancillary_Path : pathlib.Path object
        Directory containing ancillary data file(s).
    package_top_Path : pathlib.Path object
        Top-level directory of this package.
    """

    # Check for externally-provided archival version ID, construct file-search
    #  string accordingly, then determine input file(s):
    pattern = re.compile(r"\d\d-.*")
    if pattern.match(input_fpath_or_fpathpfx):
        ext_archival_versionID = input_fpath_or_fpathpfx[0:2]
        in_fp_or_fppfx = input_fpath_or_fpathpfx[3:]
    else:
        ext_archival_versionID = None
        in_fp_or_fppfx = input_fpath_or_fpathpfx
    if '*' in in_fp_or_fppfx:  # This seems to be a fpathpfx
        search_str = in_fp_or_fppfx+'*'
    elif in_fp_or_fppfx.split('.')[-1] in ["nc", "bin", "csv"]:
        search_str = in_fp_or_fppfx
    else:  # This seems to be a fpathpfx
        search_str = in_fp_or_fppfx+'*'

    input_fpath_l = glob.glob(search_str)
    n_input_fpaths = len(input_fpath_l)
    if n_input_fpaths < 1:
        raise FileNotFoundError("No input file(s) for the product generator "
                          "were found for search string {}".format(search_str))

    input_file_type, processing_plan, input_file_info = (
                                     input_file_processing_plan(input_fpath_l))

    if "concatenate" in processing_plan:
        # Concatenate the input files into a new singular (output) file:
        input_fpath = concatenate_along_dim(input_fpath_l, output_Path)
    else:
        # If not concatenating, only a single input file is supported:
        input_fpath = input_fpath_l[0]

    if "finalize" in processing_plan:
        output_fPath, granule_subset = finalize_product_file(input_fpath,
                  input_file_type, input_file_info, output_Path, output_options)
        product_fpath = str(output_fPath)
    else:
        granule_subset = False
        product_fpath = input_fpath

    if "extract_metrics" in processing_plan:
        extract_metrics(product_fpath, str(output_Path))

    Wvw_coG_fpath_l = []
    if "PRD_GEN_BASIC_ONLY" not in os.environ:
        if "create_NASA_EDS_repr_img" in processing_plan:
            import PREFIRE_PRD_GEN.product_plot_EDS_repr as EDS_repr_plot

            EDS_repr_plot.main(product_fpath, str(output_Path), granule_subset,
                               ancillary_Path)

        if "create_NASA_Worldview_COG" in processing_plan:
            import PREFIRE_PRD_GEN.product_plot_Worldview_repr as Worldview_coG

            Wvw_coG_fpath_l = Worldview_coG.main(product_fpath,
                                                 str(output_Path))

        if "create_DAAC_metadata" in processing_plan:
            if len(Wvw_coG_fpath_l) > 0:
                for Wvw_coG_fpath in Wvw_coG_fpath_l:
                    aux_input_file_info = {"assoc_nc_prd_fpath": product_fpath}
                    produce_DAAC_metadata(output_Path, Wvw_coG_fpath,
                                  input_file_type, input_file_info,
                                  ext_archival_versionID,
                                  aux_input_file_info=aux_input_file_info)

            produce_DAAC_metadata(output_Path, product_fpath,
                                  input_file_type, input_file_info,
                                  ext_archival_versionID)
