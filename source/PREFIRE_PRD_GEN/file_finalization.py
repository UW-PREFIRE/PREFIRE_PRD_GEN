"""
This program requires python version 3.6 or later, and is importable as a
python module.
"""

  # From the Python standard library:
from pathlib import Path
import datetime
import os
import json
import copy

  # From other external Python packages:
import netCDF4
import numpy as np

  # Custom utilities:
from PREFIRE_PRD_GEN.file_creation import write_data_fromspec
from PREFIRE_PRD_GEN.file_extraction import read_data_fromspec
import PREFIRE_PRD_GEN.filepaths as filepaths


superceded_global_atts = ("file_name", "granule_type", "netCDF_lib_version",
                          "UTC_of_file_creation")


def subset_data_by_atrack_mask(in_dat, atrack_mask, filespecs_json_fpath,
                               use_shared_geometry_filespecs=True):
    """
    Given a boolean mask (dimensioned 'atrack'), subset all relevant
     (atrack-dimensioned) fields in in_dat by replacing all unmasked data field
     values with fill_value).

    Input Parameters
    ----------
    in_dat : dict
        A nested Python dictionary with data groups and fields
    atrack_mask : boolean ndarray
        A boolean 1-D mask (dimensioned 'atrack'), where True preserves the
         in_dat field's value, and False sets the dat output field's value to
         fill_value.
    filespecs_json_fpath : str
        Filepath of a JSON-format file with the data specification. It is
        assumed that the keys in 'in_dat' correspond to this spec file; any
        mismatched keys in 'in_dat' are simply ignored during file creation.
    use_shared_geometry_filespecs : boolean
        (optional) If False, no shared geometry group filespec info will be
                   added.  If True (the default), add the shared geometry group
                   filespec info (stored in "shared_geometry_filespecs.json").
    """

    # Load the given product file data specification(s):
    with open(filespecs_json_fpath, 'r') as f:
        product_filespecs = json.load(f)
    full_filespecs = product_filespecs.copy()  # Init full filespec dict

    # Load any shared geometry group file data specification:
    if use_shared_geometry_filespecs:
        with open(filepaths._shared_geometry_filespec_fpath, 'r') as f:
            geometry_filespecs = json.load(f)
        full_filespecs.update(geometry_filespecs)  # Augment full filespec dict

    dat = copy.deepcopy(in_dat)

    for gn in dat:
        if ("Attributes" in gn) or ("Geometry" in gn):
            continue
        for vn in full_filespecs[gn]:
            ff_var = full_filespecs[gn][vn]
            dimn_l = ff_var["C_dimensions"]
            if ("atrack" in dimn_l) and ("bitflags" not in vn):
                shp = dat[gn][vn].shape

                n_dims = len(dimn_l)
                if n_dims == 1:
                    msk = atrack_mask
                elif n_dims == 2:
                    msk = np.broadcast_to(atrack_mask[:,np.newaxis], shp)
                elif n_dims == 3:
                    msk = np.broadcast_to(
                                      atrack_mask[:,np.newaxis,np.newaxis], shp)

                dat[gn][vn].mask = ~msk

    return dat


#--------------------------------------------------------------------------
def finalize_product_file(input_fpath, input_file_type, input_file_info,
                          output_Path, output_options):
    """Finalize a NetCDF-format product file."""

    if output_options is not None:
        outp_opts = output_options.split('|')

    product_moniker = input_file_type["product_moniker"]
    
    dat = {}

    # Read all relevant fields from the input file:
    with netCDF4.Dataset(input_fpath, 'r') as nc_ds:
        for pt_group in input_file_info["product_type_l"]:
            tmpdat = read_data_fromspec(nc_ds,
                             input_file_info["filespecs_json_fpath"], pt_group,
                        use_shared_geometry=input_file_info["use_shared_geom"],
                                 global_atts_to_ignore=superceded_global_atts,
                                 verbose=False)
            dat.update(tmpdat)  # Merge new info into this dictionary

    pass_d = {"_desc_pass": -1, "_asc_pass": 1, "_all_passes": None}

    # Process output options:
    granule_subset = False

    # For now, hardwire latitude bounds of L2+ finalized output:
    outp_opts = output_options
    if outp_opts is None:
        if product_moniker[0:2] in ["2B", "3-"]:
            outp_opts = "lat_bounds_for_all_passes: -90. to -60., 60. to 90."

    if outp_opts is not None:
        outp_opts_p = outp_opts.split('|')

        tmp = [x.split(':') for x in outp_opts_p if "lat_bounds" in x]
        if (len(tmp) > 0) and ("Geometry" in dat):
          # Subset in latitude and/or pass type (only valid for products that
          #  contain the 'Geometry' data group):
            granule_subset = True

            tokens = tmp[0][0].split("for")
            select_pass_type = None
            if len(tokens) == 2:
                select_pass_type = pass_d[tokens[1].strip()]

            lat_bounds_parts = [x.strip().split("to") for x
                                                in tmp[0][1].split(',')]
            lat_bounds = [(float(xb), float(xe)) for xb, xe in lat_bounds_parts]

            gd = dat["Geometry"]
            lat = gd["latitude"]
            atrack_dim = (lat.shape)[0]

            if select_pass_type is not None:
                sat_pass_msk = (gd["satellite_pass_type"] == select_pass_type)
            else:
                sat_pass_msk = np.full((atrack_dim,), True)

            atrack_lat_msk = np.full((atrack_dim,), False)
            for bounds in lat_bounds:
                # If any scene/FOV center in a given frame is within the latitude
                #  bounds, keep the entire frame
                lat_bnds_msk = ((lat > bounds[0]) & (lat < bounds[1]))
                atrack_lat_msk = ((np.any(lat_bnds_msk, axis=1) &
                                   sat_pass_msk) | atrack_lat_msk)  # {atrack}

            dat = subset_data_by_atrack_mask(dat, atrack_lat_msk,
                                       filepaths._final_products_filespec_fpath,
               use_shared_geometry_filespecs=input_file_info["use_shared_geom"])

    dat_gatts = dat["Global_Attributes"]
    
    # Create output filepath:
    
    date_and_time = dat_gatts["UTC_coverage_start"].split('T')
    tmp_sdate = date_and_time[0].replace('-', '')
    tmp_stime = date_and_time[1].split('.')[0].replace(':', '')
    UTC_coverage_start_intrep = tmp_sdate+tmp_stime

    if product_moniker[0:2] == "3-":
        date_and_time = dat_gatts["UTC_coverage_end"].split('T')
        tmp_sdate = date_and_time[0].replace('-', '')
        tmp_stime = date_and_time[1].split('.')[0].replace(':', '')
        UTC_coverage_end_intrep = tmp_sdate+tmp_stime

        output_fname = "PREFIRE_SAT{}_{}_{}_{}_{}.nc".format(
                        dat_gatts["spacecraft_ID"][8:], product_moniker,
                        dat_gatts["full_versionID"], UTC_coverage_start_intrep,
                        UTC_coverage_end_intrep)
    elif product_moniker[0:2] == "0-":
        output_fname = os.path.basename(input_fpath)[4:]
    else:  # 1A-, 1B-, 2B-, AUX-, ANC-
        output_fname = "PREFIRE_SAT{}_{}_{}_{}_{}.nc".format(
                        dat_gatts["spacecraft_ID"][8:], product_moniker,
                        dat_gatts["full_versionID"], UTC_coverage_start_intrep,
                        dat_gatts["granule_ID"])
    output_fPath = output_Path / output_fname

    # Set or supercede some global attribute values:

    dat["Global_Attributes"]["file_name"] = output_fname
    dat["Global_Attributes"]["granule_type"] = product_moniker
    dat["Global_Attributes"]["netCDF_lib_version"] = (
                                            netCDF4.getlibversion().split()[0])
    with open(filepaths.scipkg_prdgitv_fpath, 'r') as in_f:
        line_parts = in_f.readline().split('(', maxsplit=1)
        dat["Global_Attributes"]["provenance"] = "{}, finalization ( {}".format(
                  dat["Global_Attributes"]["provenance"], line_parts[1].strip())

    now_UTC_DT = datetime.datetime.now(datetime.timezone.utc)
    now_UTC_str = now_UTC_DT.strftime("%Y-%m-%dT%H:%M:%S.%f")
    dat["Global_Attributes"]["UTC_of_file_creation"] = now_UTC_str

    # Write finalized product file:
    write_data_fromspec(dat, str(output_fPath),
                            input_file_info["filespecs_json_fpath"],
              use_shared_geometry_filespecs=input_file_info["use_shared_geom"],
                            verbose=False)

    return (output_fPath, granule_subset)
