"""
Routines for writing a NetCDF-format file based on data file specification(s).

This program requires python version 3.6 or later, and is importable as a
python module.
"""

  # From the Python standard library:
import json
import warnings
import collections
import copy

  # From other external Python packages:
import netCDF4
import numpy as np

  # Custom utilities:
import PREFIRE_PRD_GEN.filepaths as filepaths


def _validate_shapes(dat, dat_descrs, verbose=False):
    """ Helper routine to validate the shapes. 'dat' is a nested dictionary
    containing the data arrays; 'dat_descrs' is the nested dictionary loaded
    from the JSON-format data specification file.

    This routine makes sure the arrays in 'dat' have consistent shapes, as
    indicated by the dimension names in 'dat_descrs'.

    If 'verbose', it prints a summary of the dimension names it found.

    Returns a dictionary, with the keys equal to the dimension names
    (strings) and the values equal to the dimension lengths (integers)
    """

    dim_shapes = collections.defaultdict(set)

    for g_name in dat:

        if g_name not in dat_descrs:
            if verbose:
                warnings.warn(
                        "Group {} is not present in dat_descrs, omitting it". \
                                                                format(g_name))
            continue

        if "_Attributes" not in g_name:
            for v_name in dat[g_name]:
                if v_name not in dat_descrs[g_name]:
                    if verbose:
                        warnings.warn("Var {} is not present in dat_descrs, "
                                      "omitting it".format(v_name))
                    continue
                for d, dim in enumerate(
                                   dat_descrs[g_name][v_name]["C_dimensions"]):
                    if "!U!" in dim:  # Unlimited dimension
                        actual_dim = dim.replace("!U!", '')
                        dim_shapes[actual_dim].add(0)
                    elif "nchars" in dim:  # Character-string dimension
                        # Only works for 1-D arrays of "strings"
                        dim_shapes[dim].add(len(dat[g_name][v_name][0]))
                    else:  # Regular dimension
                        dim_shapes[dim].add(dat[g_name][v_name].shape[d])

    for dim in dim_shapes:
        if len(dim_shapes[dim]) > 1:
            raise ValueError("Dimension {} size is not consistent: {} ". \
                             format(dim, str(dim_shapes[dim])))
        if verbose:
            print("Dimension {}, length {}".format(dim, str(dim_shapes[dim])))

    # Convert the sets to scalar integers:
    out_dim_shapes = {}
    for dim in dim_shapes:
        out_dim_shapes[dim] = dim_shapes[dim].pop()

    return out_dim_shapes


def _write_generic(nc, dat, dat_descrs, verbose=False,
                   omit_missing_vars=False):
    """
    Generic NetCDF-format file writer.

    Parameters
    ----------
    nc : netCDF4.Dataset
        the netCDF4.Dataset object, referencing the desired file
    dat : dict
        a nested Python dictionary with data to write
    dat_descrs : dict
        the loaded JSON data description information.
    omit_missing_vars : bool, optional
        whether to omit variables from the output file that are present in the
        file specs but missing in the output data dictionary.
        The default is False.
    """

    # Method:
    # 1) Validate and identify dimension lengths of variables within 'dat',
    #     using a helper function.
    # 2) Find the actual groups defined in 'dat'.
    # 3) For each group defined in 'dat', copy all variables and
    #     global/group attributes from 'dat_descrs'.
    # 4) For each variable and global/group attribute, get the data inside
    #     'dat'.
    # 5) Then proceed to define and populate the NetCDF file:
    #      Set any global attributes
    #      Create the first group
    #        Set any group attributes
    #        Create each variable
    #          Set any dimensions, if needed.
    #          Write the data, if present (otherwise, remains FillValue)

    v_atts_to_ignore = ("np_dtype", "C_dimensions", "fill_value", "compression",
                        "fullname", "v_name", "g_name", "data", "F_dimensions")
    v_atts_same_dtype_as_var = ["flag_values", "valid_min", "valid_max"]

    dim_shapes = _validate_shapes(dat, dat_descrs, verbose=verbose)

    groups_to_write = set()
    g_name = "Global_Attributes"
    if g_name in dat_descrs:
        groups_to_write.add(g_name)
    for g_name in dat:
        if g_name in dat_descrs:
            groups_to_write.add(g_name)
            ga_name = g_name+"_Group_Attributes"
            if ga_name in dat_descrs:
                groups_to_write.add(ga_name)

    vars_to_write = []
    gatts_to_write = []
    gpatts_to_write = []
    for g_name in groups_to_write:
        if g_name == "Global_Attributes":  # Contains global attribute(s)
            for gatt_name in dat_descrs[g_name]:
                new_gatt = {}
                # create a new dictionary so we can add the data key to it
                new_gatt.update(dat_descrs[g_name][gatt_name])
                new_gatt["fullname"] = gatt_name
                new_gatt["g_name"] = g_name
                new_gatt["gatt_name"] = gatt_name
                if dat_descrs[g_name][gatt_name]["value"] == "!to_be_set!":
                    a_val = None
                else:
                    a_val = dat_descrs[g_name][gatt_name]["value"]
                if g_name in dat:
                    if gatt_name in dat[g_name]:
                        a_val = dat[g_name][gatt_name]
                if a_val is not None:
                    new_gatt["data"] = a_val
                else:
                    raise ValueError("the value of global attribute {} "
                                     "has not been set.".format(gatt_name))
                gatts_to_write.append(new_gatt)
        elif "_Group_Attributes" in g_name:  # Contains group attribute(s)
            for gpatt_name in dat_descrs[g_name]:
                actual_gname = g_name.replace("_Group_Attributes", '')
                new_gpatt = {}
                # create a new dictionary so we can add the data key to it
                new_gpatt.update(dat_descrs[g_name][gpatt_name])
                new_gpatt["fullname"] = gpatt_name
                new_gpatt["g_name"] = actual_gname
                new_gpatt["gpatt_name"] = gpatt_name
                if dat_descrs[g_name][gpatt_name]["value"] == "!to_be_set!":
                    a_val = None
                else:
                    a_val = dat_descrs[g_name][gpatt_name]["value"]
                if g_name in dat:
                    if gpatt_name in dat[g_name]:
                        a_val = dat[g_name][gpatt_name]
                if a_val is not None:
                    new_gpatt["data"] = a_val
                else:
                    raise ValueError("the value of group ({}) attribute {} "
                                     "has not been set.".format(actual_gname,
                                                                gpatt_name))
                gpatts_to_write.append(new_gpatt)
        else:  # This is a regular group
            for v_name in dat_descrs[g_name]:
                if (v_name not in dat[g_name]) and omit_missing_vars:
                    continue

                new_var = {}
                  # create a new dictionary so we can add the data key to it
                new_var.update(dat_descrs[g_name][v_name])
                new_var["fullname"] = g_name+'/'+v_name
                new_var["g_name"] = g_name
                new_var["v_name"] = v_name

                  # Ensure the dtype of these match the var dtype:
                for aname in v_atts_same_dtype_as_var:
                    if aname in new_var:
                        try:
                            n = len(new_var[aname])
                            new_var[aname] = np.array(new_var[aname],
                                                      dtype=new_var["np_dtype"])
                        except TypeError:  # Probably a scalar
                            tmp = np.array([new_var[aname]],
                                                      dtype=new_var["np_dtype"])
                            new_var[aname] = tmp[0]

                if v_name in dat[g_name]:
                    new_var["data"] = dat[g_name][v_name]
                vars_to_write.append(new_var)

    for gatt in gatts_to_write:
        nc.setncattr(gatt["fullname"], gatt["data"])

    g_inst = {}  # Init dict of Group class instances

    for gpatt in gpatts_to_write:
        if gpatt["g_name"] not in nc.groups:
            g_inst[gpatt["g_name"]] = nc.createGroup(gpatt["g_name"])
        g_inst[gpatt["g_name"]].setncattr(gpatt["fullname"], gpatt["data"])

    for var in vars_to_write:

        if var["g_name"] not in nc.groups:
            g_inst[var["g_name"]] = nc.createGroup(var["g_name"])
    
        dims_exist = True
        for d, dimen in enumerate(var["C_dimensions"]):
            if dimen not in nc.dimensions:
                # Create the dimension, if it is in dim_shapes (which is what
                #  was able to be defined from the actual input data):
                if dimen in dim_shapes:
                    nc.createDimension(dimen, dim_shapes[dimen])
                else:
                    warnings.warn("cannot write variable {} to file, as no "
                                  "data for the variable was present to "
                                  "determine the size of at least one of its "
                                  "dimensions".format(var["fullname"]))
                    dims_exist = False
                    break

        if dims_exist:
            kw_d0 = {"dimensions": var["C_dimensions"]}

            if "fill_value" in var:
                kw_d0["fill_value"] = var["fill_value"]

            if "compression" in var:
                if len(var["compression"]) > 2 and "float" in var["np_dtype"]:
                    compression_level, use_shuffle, lossy = var["compression"]
                    if "ndecplc" in lossy:  # e.g., "ndecplc_3"
                        # Use 'least_significant_digit' keyword argument, which
                        #  specifies the power of ten of the smallest decimal
                        #  place in the data that is a reliable value
                        # (performed by Python code, not the underlying NetCDF
                        #  library)
                        kw_d0["least_significant_digit"] = int(lossy[8:])
                    elif "nsigdig" in lossy:  # e.g., "nsigdig_6"
                        # Use 'significant_digits' keyword argument, which
                        #  specifies the number of significant digits
                        #  independent of the magnitude of the variable
                        # (requires NetCDF C library version >= 4.9.0)
                        kw_d0["significant_digits"] = int(lossy[8:])
                else:
                    compression_level, use_shuffle = var["compression"]
#                kw_d0["quantize_mode"] = "GranularBitRound"

                kw_d0["complevel"] = compression_level
                kw_d0["shuffle"] = use_shuffle

                kw_newer_d = copy.deepcopy(kw_d0)
                kw_newer_d["compression"] = "zlib"

                kw_old_d = copy.deepcopy(kw_d0)
                kw_old_d["zlib"] = True

                try:
                    v_obj = nc.createVariable(var["fullname"],
                                              var["np_dtype"], **kw_newer_d)
                except:  # Try a method that may work for older software
                    v_obj = nc.createVariable(var["fullname"],
                                              var["np_dtype"], **kw_old_d)

                if ("fill_value" not in var) and (var["np_dtype"] == "S1"):
                    # Only works for 1-D arrays of "strings"
                    tmp_data = copy.deepcopy(var["data"])
                    var["data"] = np.array([list(x) for x in tmp_data])

            for aname in var:
                if aname not in v_atts_to_ignore:
                    v_obj.setncattr(aname, var[aname])

            if "data" in var:
                nc[var["fullname"]][...] = var["data"]


def write_data_fromspec(dat, output_fpath, filespecs_json_fpath,
                        use_shared_geometry_filespecs=True, verbose=False,
                        omit_missing_vars=False):
    """
    Write data to file using given specification(s).

    Input Parameters
    ----------
    dat : dict
        A nested Python dictionary with data groups and fields
    output_fpath : str
        Filepath of the output file that will be created
    filespecs_json_fpath : str
        Filepath of a JSON-format file with the data specification. It is
        assumed that the keys in 'dat' correspond to this spec file; any
        mismatched keys in 'dat' are simply ignored during file creation.
    use_shared_geometry_filespecs : boolean
        (optional) If False, no shared geometry group filespec info will be
                   added.  If True (the default), add the shared geometry group
                   filespec info (stored in "shared_geometry_filespecs.json").
    verbose : bool
        (optional) Should verbose details be sent to stdout/stderr?
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

    with netCDF4.Dataset(output_fpath, 'w', clobber=True) as nc_ds:
        _write_generic(nc_ds, dat, full_filespecs, verbose=verbose,
                       omit_missing_vars=omit_missing_vars)
