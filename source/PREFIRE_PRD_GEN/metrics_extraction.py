"""
This program requires python version 3.6 or later, and is importable as a
Python module.
"""

  # From the Python standard library:
import os
import json

  # From other external Python packages:
import netCDF4
import numpy as np

  # Custom utilities:
from PREFIRE_PRD_GEN.file_creation import write_data_fromspec
from PREFIRE_PRD_GEN.file_read import load_all_nc4_global_atts, \
                                      get_PREFIRE_Lx_field
import PREFIRE_PRD_GEN.filepaths as filepaths


cfg_d = {}

cfg_d["0-PAYLOAD-TLM"] = {
        "spec_group": "TIRS_L0_Data",
        "metrics": ["encoder_stats"],
        "encoder_stats": ["stats", "encoder_pos", "preserve_no_dims"]}

cfg_d["0-BUS-TLM"] = {
        "spec_group": "SCbus_L0_Data",
        "metrics": ["granule_edge_ctime"],
        "granule_edge_ctime": ["copy", "granule_edge_ctime",
                               "preserve_all_dims"]}

cfg_d["0-ORB-REC"] = {
        "spec_group": "Orbit_L0_Data",
        "metrics": ["granule_edge_ctime"],
        "granule_edge_ctime": ["copy", "granule_edge_ctime",
                               "preserve_all_dims"]}

cfg_d["1A-RAD"] = {
        "spec_group": "NonGeoLoc_Radiance",
        "metrics": ["n_nonnom_radobsq"],
        "n_nonnom_radobsq": ["count_gt_this_value", "radiance_quality_flag",
                             "preserve_atrack", 0]}

cfg_d["1B-RAD"] = {
        "spec_group": "Radiance",
        "metrics": ["n_nonnom_geoloc"],
        "n_nonnom_geoloc": ["count_gt_this_value", "geoloc_quality_bitflags",
                            "preserve_atrack", 0, "Geometry"]}

cfg_d["AUX-MET"] = {
        "spec_group": "Aux-Met",
        "metrics": ["n_missing_skin_temp"],
        "n_missing_skin_temp": ["count_missing", "skin_temp",
                                "preserve_atrack"]}

cfg_d["AUX-SAT"] = {
        "spec_group": "Aux-Sat",
        "metrics": ["n_merged_seaice_f_dsrc", "n_merged_snow_f_dsrc"],
        "n_merged_seaice_f_dsrc": ["count_val_to_dim_entry",
                                   "merged_seaice_final_data_source",
                                   "preserve_atrack"],
        "n_merged_snow_f_dsrc": ["count_val_to_dim_entry",
                                 "merged_snow_final_data_source",
                                 "preserve_atrack"]}

cfg_d["2B-MSK"] = {
        "spec_group": "Msk",
        "metrics": ["n_missing_cloud_mask"],
        "n_missing_cloud_mask": ["count_missing", "cloud_mask",
                                 "preserve_atrack"]}

cfg_d["2B-SFC"] = {
        "spec_group": "Sfc",
        "metrics": ["n_missing_sfc_sp_emis"],
        "n_missing_sfc_sp_emis": ["count_missing", "sfc_spectral_emis",
                                  "preserve_atrack"]}

cfg_d["2B-FLX"] = {
        "spec_group": "Flx",
        "metrics": ["n_missing_olr", "n_missing_spectral_flux"],
        "n_missing_olr": ["count_missing", "olr", "preserve_atrack"],
        "n_missing_spectral_flux": ["count_missing", "spectral_flux",
                                    "preserve_atrack"]}

cfg_d["2B-ATM"] = {
        "spec_group": "Atm",
        "metrics": ["n_missing_cwv", "n_missing_T_profile",
                    "n_missing_wv_profile", "n_missing_surface_T"],
        "n_missing_cwv": ["count_missing", "cwv", "preserve_atrack"],
        "n_missing_T_profile": ["count_missing", "T_profile",
                                "preserve_atrack"],
        "n_missing_wv_profile": ["count_missing", "wv_profile",
                                 "preserve_atrack"],
        "n_missing_surface_T": ["count_missing", "surface_T",
                                "preserve_atrack"]}

cfg_d["2B-CLD"] = {
        "spec_group": "Cld",
        "metrics": ["n_missing_cloudtop_p", "n_missing_cloud_tau",
                    "n_missing_cloud_d_eff"],
        "n_missing_cloudtop_p": ["count_missing", "cloudtop_pressure",
                                 "preserve_atrack"],
        "n_missing_cloud_tau": ["count_missing", "cloud_tau",
                                "preserve_atrack"],
        "n_missing_cloud_d_eff": ["count_missing", "cloud_d_eff",
                                  "preserve_atrack"]}

cfg_d["3-SFC-SORTED-ALLSKY"] = {
        "spec_group": "Sfc-Sorted",
        "metrics": ["n_missing_emis_mean"],
        "n_missing_emis_mean": ["count_missing", "emis_mean",
                                "preserve_latitude"]}


#------------------------------------------------------------------------------
def granule_info(p_fpath):
    """Determines top-level info about the given granule."""

    try:
        with netCDF4.Dataset(p_fpath, 'r') as nc_ds:
            try:
                sat_num = int(nc_ds.sensor_ID[5])
            except:
                sat_num = int(nc_ds.spacecraft_ID[8])
            granule_IDstr = nc_ds.granule_ID

            granule_type = nc_ds.granule_type
            if granule_type == "1B-NLRAD":
                granule_type = "1B-RAD"
    except:
        msg = f"Input file ({p_fpath}) not found or unreadable"
        raise RuntimeError(msg)

    info_d = {"sat_num": sat_num, "granule_IDstr": granule_IDstr,
              "granule_type": granule_type}

    return info_d


#------------------------------------------------------------------------------
def extract_metrics(p_fpath, outp_path):
    """Main driver for extraction of metrics."""

    # Input parameters:
    #   p_fpath  :  Filepath of product datafile
    #   outp_path  :  path (preferably an absolute one) to write the output
    #                 file(s) to

    dat = {}  # Initialize output data dictionary

    # Determine granule info:
    info_d = granule_info(p_fpath)

    nc_ds = netCDF4.Dataset(p_fpath, 'r')

    granule_type = info_d["granule_type"]

    # Extract all global attributes (most/all of these will simply be copied
    #  into the output metrics file), and initialize global_atts dictionary:
    global_atts = load_all_nc4_global_atts(nc_ds)

    global_atts["summary"] = ("Monitoring metrics for the PREFIRE {} "
                              "product.".format(granule_type))

    # Load the metrics file data specification(s):
    filespecs_json_fpath = filepaths._extracted_metrics_filespec_fpath
    with open(filespecs_json_fpath, 'r') as f:
        metrics_filespecs = json.load(f)

    mcfg_d = cfg_d[granule_type]
    spgrp = mcfg_d["spec_group"]

    dat[spgrp] = {}  # Initialize ouput group

    for metric in mcfg_d["metrics"]:
        outp_dtype = metrics_filespecs[spgrp][metric]["np_dtype"]

        if mcfg_d[metric][2].split('_')[1] == "atrack":
            along_track_np_index_range = ("atrack", 0,
                                               nc_ds.dimensions["atrack"].size)
        else:
            along_track_np_index_range = None

        if mcfg_d[metric][0] == "count_missing":
            # Extract the field:
            field_name = mcfg_d[metric][1]
            field = get_PREFIRE_Lx_field(nc_ds, spgrp, field_name,
                                         along_track_np_index_range)

            # Count along a "dimension" to keep:
            try:
                dimname_to_keep = mcfg_d[metric][2].split('_')[1]
                ii = (nc_ds.groups[spgrp].variables[field_name].dimensions.
                      index(dimname_to_keep))
                ax_t = tuple(x for x in range(len(field.shape)) if x != ii)
            except:
                ii = None
            count = np.count_nonzero(np.ma.getmaskarray(field), axis=ax_t)
            dat[spgrp][metric] = count.astype(outp_dtype)

        elif mcfg_d[metric][0] == "stats":
            # Extract the field:
            field_name = mcfg_d[metric][1]
            field = get_PREFIRE_Lx_field(nc_ds, spgrp, field_name,
                                         along_track_np_index_range)

            # (0)minimum, (1)mean, (2)median, (3)maximum, (4)standard deviation
            field_min = field.min()
            field_mean = field.mean()
            field_median = np.median(field)
            field_max = field.max()
            field_std = np.std(field)
            dat[spgrp][metric] = np.array([field_min, field_mean, field_median,
                                    field_max, field_std])

        elif mcfg_d[metric][0] == "copy":
            # Extract the field:
            field_name = mcfg_d[metric][1]
            dat[spgrp][metric] = get_PREFIRE_Lx_field(nc_ds, spgrp, field_name,
                                               along_track_np_index_range)

        elif mcfg_d[metric][0] == "count_gt_this_value":
            # Extract the field:
            field_name = mcfg_d[metric][1]
            fldgrp = spgrp
            if len(mcfg_d[metric]) > 4:  # Another group name is specified
                fldgrp = mcfg_d[metric][4]
            field = get_PREFIRE_Lx_field(nc_ds, fldgrp, field_name,
                                         along_track_np_index_range)

            # Count along a "dimension" to keep:
            try:
                dimname_to_keep = mcfg_d[metric][2].split('_')[1]
                ii = (nc_ds.groups[fldgrp].variables[field_name].dimensions.
                      index(dimname_to_keep))
                ax_t = tuple(x for x in range(len(field.shape)) if x != ii)
            except:
                ii = None
            count = np.count_nonzero(field > mcfg_d[metric][3], axis=ax_t)
            dat[spgrp][metric] = count.astype(outp_dtype)

        elif mcfg_d[metric][0] == "count_val_to_dim_entry":
            # Extract the field:
            field_name = mcfg_d[metric][1]
            field = get_PREFIRE_Lx_field(nc_ds, spgrp, field_name,
                                         along_track_np_index_range)

            maxv = field.max()
            if maxv > 30:
                raise RuntimeError("The maximum value of this field is "
                                   "unsuited for this operation.")

            # Count along a "dimension" to keep:
            try:
                dimname_to_keep = mcfg_d[metric][2].split('_')[1]
                ii = (nc_ds.groups[spgrp].variables[field_name].dimensions.
                      index(dimname_to_keep))
                ax_t = tuple(x for x in range(len(field.shape)) if x != ii)
                size_dtk = nc_ds.dimensions[dimname_to_keep].size
                count = np.zeros((int(maxv), size_dtk), dtype=outp_dtype)
            except:
                ii = None

            for i in range(field.max()):
                count[0,:] = np.count_nonzero(field == i, axis=ax_t)
            dat[spgrp][metric] = count

      # Add global attributes to output dictionary:
    dat["Global_Attributes"] = global_atts
                                 
    fn = os.path.basename(p_fpath)
    output_fpath = os.path.join(outp_path, fn[:-2]+"metrics.nc")

    write_data_fromspec(dat, output_fpath, filespecs_json_fpath,
                        use_shared_geometry_filespecs=False, verbose=False)

    nc_ds.close()

