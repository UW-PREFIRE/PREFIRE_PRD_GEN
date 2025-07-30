"""
Routines for extracting information from a NetCDF-format file based on data
 file specification(s).

This program requires python version 3.6 or later, and is importable as a
python module.
"""

  # From the Python standard library:
import sys
import json
import collections

  # From other external Python packages:

  # Custom utilities:
import PREFIRE_PRD_GEN.filepaths as filepaths
from PREFIRE_tools.utils.numeric import contains_NaN_or_Inf, replace_NaN_or_Inf


def _dict_from_relevant_data(nc_ds, full_filespecs, product_type,
                             use_shared_geometry, global_atts_to_ignore,
                             verbose):
    """Populate a dictionary with relevant data from 'nc_ds'."""

    dat = {}

    inpf_group_names = [x.name for x in nc_ds.groups.values()]
    outp_group_names = [x for x in full_filespecs.keys() \
                                                     if x != "_JSON_COMMENTS_"]

    # Process any global attributes:
    global_att_names = nc_ds.ncattrs()
    if len(global_att_names) >= 1:
        g_name = "Global_Attributes"
        dat[g_name] = {}
        for a_name in full_filespecs[g_name]:
             if (full_filespecs[g_name][a_name]["value"] == "!to_be_set!") and \
                     (a_name not in global_atts_to_ignore):
                 try:   # REMOVE this line later
                     tmp_o = nc_ds.getncattr(a_name)
                     if contains_NaN_or_Inf(tmp_o):
                         dat[g_name][a_name] = replace_NaN_or_Inf(tmp_o,
                                                            warn_with_moniker=
                                           f"{{group {g_name}, att {a_name}}}")
                     else:
                         dat[g_name][a_name] = tmp_o
                 except:   # REMOVE this line (and the two below) later
                     if a_name in ["orbit_sim_version", "SRF_NEdR_version"]:
                         dat[g_name][a_name] = "unknown"

    # Determine group names relevant to this product type:
    relevant_group_names = [product_type]
    ptga = product_type+"_Group_Attributes"
    if ptga in outp_group_names:
        relevant_group_names.append(ptga)
    if use_shared_geometry:
        if "NonGeoLoc_" in product_type or product_type == "Cal_Artifacts":
            ptg = "ProtoGeometry"
        else:
            ptg = "Geometry"
        relevant_group_names.append(ptg)
        ptga = ptg+"_Group_Attributes"
        if ptga in outp_group_names:
            relevant_group_names.append(ptga)

    # Process relevant groups:
    for g_name in relevant_group_names:
        gname_has_ncdata = (g_name in outp_group_names)
        if gname_has_ncdata:
            dat[g_name] = {}
        else:
            continue  # Nothing more to do this iteration

        if "Group_Attributes" in g_name:  # group attributes
            for a_name in full_filespecs[g_name]:
                if full_filespecs[g_name][a_name]["value"] == "!to_be_set!":
                    real_gname = g_name.replace("_Group_Attributes", '')
                    tmp_o = nc_ds.groups[real_gname].getncattr(a_name)
                    if contains_NaN_or_Inf(tmp_o):
                        dat[g_name][a_name] = replace_NaN_or_Inf(tmp_o,
                                                            warn_with_moniker=
                                           f"{{group {g_name}, att {a_name}}}")
                    else:
                        dat[g_name][a_name] = tmp_o

        else:  # group variables
            for v_name in full_filespecs[g_name]:
                if "fill_value" in full_filespecs[g_name][v_name]:
                    fill_value = full_filespecs[g_name][v_name]["fill_value"]
                else:
                    fill_value = None
                if "optional" in full_filespecs[g_name][v_name]:
                    if full_filespecs[g_name][v_name]["optional"] == "yes":
                        try:
                            tmp_o = nc_ds.groups[g_name].variables[v_name][...]
                            if contains_NaN_or_Inf(tmp_o):
                                dat[g_name][v_name] = replace_NaN_or_Inf(tmp_o,
                                      fill_value=fill_value, warn_with_moniker=
                                           f"{{group {g_name}, var {v_name}}}")
                            else:
                                dat[g_name][v_name] = tmp_o
                        except:
                            pass
                else:
                    tmp_o = nc_ds.groups[g_name].variables[v_name][...]
                    if contains_NaN_or_Inf(tmp_o):
                        dat[g_name][v_name] = replace_NaN_or_Inf(tmp_o,
                                     fill_value=fill_value, warn_with_moniker=
                                           f"{{group {g_name}, var {v_name}}}")
                    else:
                        dat[g_name][v_name] = tmp_o
    return dat


def read_data_fromspec(nc_ds, filespecs_json_fpath, product_type,
                       use_shared_geometry=True, global_atts_to_ignore=list(),
                       verbose=False):
    """
    Read data from file using given specification(s).

    Input Parameters
    ----------
    nc_ds : netCDF4.Dataset object
        A Dataset object already associated with the NetCDF-format input file
    filespecs_json_fpath : str
        Filepath of a JSON-format file with the data specification. It is
        assumed that all (or at least some subset) of the keys/fields in
        'nc_ds' correspond to this spec file; any irrelevant mismatched
        keys in 'nc_ds' file will simply be ignored.
    product_type : str
        Group name within 'nc_ds'
    use_shared_geometry : boolean
        (optional) If False, no shared geometry group filespec info will be
                   added.  If True (the default), add the shared geometry group
                   filespec info (stored in "shared_geometry_filespecs.json").
    global_atts_to_ignore : tuple or list of strings (can be an empty sequence)
        (optional) Do not read any existing global attribute values with the
                   given names.
    verbose : bool
        (optional) Should verbose details be sent to stdout/stderr?
    """

    # Load the given product file data specification(s):
    with open(filespecs_json_fpath, 'r') as f:
        product_filespecs = json.load(f)
    full_filespecs = product_filespecs.copy()  # Init full filespec dict

    # Load any shared geometry group file data specification:
    if use_shared_geometry:
        with open(filepaths._shared_geometry_filespec_fpath, 'r') as f:
            geometry_filespecs = json.load(f)
        full_filespecs.update(geometry_filespecs)  # Augment full filespec dict

    # Populate a dictionary with relevant data from 'nc_ds':
    dat = _dict_from_relevant_data(nc_ds, full_filespecs, product_type,
                                   use_shared_geometry, global_atts_to_ignore,
                                   verbose)

    return dat
