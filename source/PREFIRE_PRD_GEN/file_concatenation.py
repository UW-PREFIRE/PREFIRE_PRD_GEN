"""
This program requires python version 3.6 or later, and is importable as a
python module.
"""

  # From the Python standard library:
from pathlib import Path
import os
import sys
import collections

  # From other external Python packages:
import numpy as np
import netCDF4

  # Custom utilities:
from PREFIRE_tools.utils.time import ctime_to_UTC_DT


def chunk_sort_key(chunk_fpath):
    """ """
    fn = os.path.basename(chunk_fpath)
    return '_'.join( fn.split('_')[-4:-2])


def concatenate_along_dim(raw_fpath_l, output_Path):
    """
    Concatenate multiple NetCDF-format files along a given named dimension.

    The dimension name and chunking information are expected to be encoded
    in the names of the input files.  A single new composite NetCDF-format
    file (with a name based on the output path and the input filepath prefix)
    results from this operation.

    Input Parameters
    ----------
    raw_fpath_l : list of str
        List of filepaths of the NetCDF-format chunkfiles to be concatenated.
    output_Path : pathlib Path object
        Directory of the output NetCDF-format file. 
    """

    if len(raw_fpath_l) < 1:
        raise ValueError("no chunk filepaths were provided")
    elif len(raw_fpath_l) == 1:
        raise ValueError("only one chunk file was provided, cannot concatenate")

    sorted_fpath_l = sorted(raw_fpath_l, key=chunk_sort_key)

    chunk_md = [tuple(((x.split('.')[-2]).split('-')[-1]).split('_')) for x in \
                sorted_fpath_l]

    if len(set([x[0] for x in chunk_md])) != 1:
        raise ValueError("the dimension name (in the filenames) is not the "
                         "same for all chunks")
    elif len(set([x[4] for x in chunk_md])) != 1:
        raise ValueError("the max dimension value (in the filenames) is not "
                         "the same for all chunks")

    concatenation_dim_name = chunk_md[0][0]

    chunk_edges = np.empty((len(chunk_md)*2), dtype=np.int32)
    i = 0
    for x in chunk_md:
        chunk_edges[i], chunk_edges[i+1] = int(x[1]), int(x[2])
        i += 2

    n_cdim = int(chunk_md[0][4].replace('f', ''))    
    if chunk_edges[0] != 0:
        raise ValueError("the initial chunk index is not zero")
    else:
        if chunk_edges[-1] != n_cdim-1:
            raise ValueError("the final chunk index is not "+str(n_cdim-1))

    if not np.all(np.diff(chunk_edges)[1::2] == 1):
        raise RuntimeError("the chunks are not all contiguous")

    # Open all dataset files:
    ds_l = []
    for sorted_fpath in sorted_fpath_l:
        ds_l.append(netCDF4.Dataset(sorted_fpath))

    ds0 = ds_l[0]  # Single out the first dataset for copying global atts, etc.

    #== Determine which items should be concatenated and which should only
    #    be copied.  Additionally, determine what dimensions need to be defined
    #    in the output file.

    copyonly_items_gl = []
    concat_items_gl = []
    ds_dim_info = collections.defaultdict(set)
    
      # For the "root" group:
    g_copy_d = {"g_name": '/'}
    g_copy_d["v_info"] = []
    g_copy_d["a_info"] = ds0.ncattrs()
    g_concat_d = {"g_name": '/'}
    g_concat_d["v_info"] = []
    for v_name in ds0.variables:
        v_dim_names = ds0.variables[v_name].dimensions
        for idim, dim_name in enumerate(v_dim_names):
            ds_dim_info[dim_name].add((ds0.variables[v_name].shape)[idim])

        varobj = ds0.variables[v_name]
        if concatenation_dim_name in v_dim_names:
            v_cdim_type = (len(v_dim_names)*10+
                                     v_dim_names.index(concatenation_dim_name))
            g_concat_d["v_info"].append((varobj, v_cdim_type))
        else:
            g_copy_d["v_info"].append(varobj)
    copyonly_items_gl.append(g_copy_d)
    concat_items_gl.append(g_concat_d)

      # For any other groups:
    for g_name in ds0.groups:
        if len(ds0.groups[g_name].groups) > 0:
            raise RuntimeError("nested groups are not yet supported for "
                               "concatenation")
        g_copy_d = {"g_name": g_name}
        g_copy_d["v_info"] = []
        g_copy_d["a_info"] = ds0.groups[g_name].ncattrs()
        g_concat_d = {"g_name": g_name}
        g_concat_d["v_info"] = []
        for v_name in ds0.groups[g_name].variables:
            v_dim_names = ds0.groups[g_name].variables[v_name].dimensions
            for idim, dim_name in enumerate(v_dim_names):
                ds_dim_info[dim_name].add(
                            (ds0.groups[g_name].variables[v_name].shape)[idim])

            # "v_info" is a list of 'netCDF4.Variable' objects:
            varobj = ds0.groups[g_name].variables[v_name]
            if concatenation_dim_name in v_dim_names:
                v_cdim_type = (len(v_dim_names)*10+
                                     v_dim_names.index(concatenation_dim_name))
                g_concat_d["v_info"].append((varobj, v_cdim_type))
            else:
                g_copy_d["v_info"].append(varobj)
        copyonly_items_gl.append(g_copy_d)
        concat_items_gl.append(g_concat_d)
    
    #== Determine the output filepath, and create/write it:

    fntoks = (os.path.basename(sorted_fpath_l[0]).split('.')[0]).split('-')
    fn = f"{fntoks[1]}-{fntoks[2]}.nc"
    output_fPath = output_Path / fn

    with netCDF4.Dataset(str(output_fPath), 'w', clobber=True) as outp_ds:

        # Create/set all dimensions:
        for dim_name, dim_vset in ds_dim_info.items():
            if dim_name != concatenation_dim_name:
                outp_ds.createDimension(dim_name, list(dim_vset)[0])
        outp_ds.createDimension(concatenation_dim_name, n_cdim)

        # Copy global attributes all at once via dictionary:
        outp_ds.setncatts(ds0.__dict__)

        ## Copy all group attributes and copy-only variables (and their
        ##  attributes):
        for item in copyonly_items_gl:
            g_name = item["g_name"]
            if g_name == '/':  # "root" group variables
                
                for varobj in item["v_info"]:
                    v_name = varobj.name

                    # Always use shuffle and "zlib" compression:
                    try:
                        outp_ds.createVariable(v_name,
                               varobj.datatype, dimensions=varobj.dimensions,
                               compression="zlib", complevel=4, shuffle=True)
                    except:
                        # Try to use the older 'zlib' argument:
                        outp_ds.createVariable(v_name,
                               varobj.datatype, dimensions=varobj.dimensions,
                               zlib=True, complevel=4, shuffle=True)

                    # Copy variable attributes all at once via dictionary
                    outp_ds[v_name].setncatts(ds0[v_name].__dict__)

                    # Copy the data values to the new variable:
                    outp_ds[v_name][...] = ds0[v_name][...]

            else:  # Other group variables and group attributes
                outp_ds.createGroup(g_name)  # Create each group here
 
                # Copy any group attributes all at once via dictionary:
                outp_ds.groups[g_name].setncatts(ds0.groups[g_name].__dict__)
                
                for varobj in item["v_info"]:
                    v_name = varobj.name

                    # Always use shuffle and "zlib" compression:
                    try:
                        outp_ds.groups[g_name].createVariable(v_name,
                               varobj.datatype, dimensions=varobj.dimensions,
                               compression="zlib", complevel=4, shuffle=True)
                    except:
                        # Try to use the older 'zlib' argument:
                        outp_ds.groups[g_name].createVariable(v_name,
                               varobj.datatype, dimensions=varobj.dimensions,
                               zlib=True, complevel=4, shuffle=True)

                    # Copy variable attributes all at once via dictionary
                    outp_ds.groups[g_name][v_name].setncatts(
                                           ds0.groups[g_name][v_name].__dict__)

                    # Copy the data values to the new variable:
                    outp_ds.groups[g_name][v_name][...] = (
                                               ds0.groups[g_name][v_name][...])

        ## Concatenate all relevant variables (copying their attributes):
        for item in concat_items_gl:
            g_name = item["g_name"]
            if g_name == '/':  # "root" group variables
                
                for varobj, _ in item["v_info"]:
                    v_name = varobj.name
                    v_dimn = varobj.dimensions
 
                    # Always use shuffle and "zlib" compression:
                    try:
                        outp_ds.createVariable(v_name,
                                 varobj.datatype, dimensions=v_dimn,
                                 compression="zlib", complevel=4, shuffle=True)
                    except:
                        # Try to use the older 'zlib' argument:
                        outp_ds.createVariable(v_name,
                                          varobj.datatype, dimensions=v_dimn,
                                          zlib=True, complevel=4, shuffle=True)

                    # Copy variable attributes all at once via dictionary
                    outp_ds[v_name].setncatts(ds0[v_name].__dict__)

            else:  # Other group variables and group attributes
                
               for varobj, _ in item["v_info"]:
                    v_name = varobj.name
                    v_dimn = varobj.dimensions

                    # Always use shuffle and "zlib" compression:
                    try:
                        outp_ds.groups[g_name].createVariable(v_name,
                                 varobj.datatype, dimensions=v_dimn,
                                 compression="zlib", complevel=4, shuffle=True)
                    except:
                        # Try to use the older 'zlib' argument:
                        outp_ds.groups[g_name].createVariable(v_name,
                                          varobj.datatype, dimensions=v_dimn,
                                          zlib=True, complevel=4, shuffle=True)

                    # Copy variable attributes all at once via dictionary
                    outp_ds.groups[g_name][v_name].setncatts(
                                           ds0.groups[g_name][v_name].__dict__)

            # Fill the data values of the new output file variables.  Loop
            #  more slowly over the datafiles (than the variables that need
            #  concatenated) in order to avoid unnecessary "drive hopping":
            ic = -2
            for ds in ds_l:
                ic += 2
                ibc = chunk_edges[ic]      # NumPy indices
                iec = chunk_edges[ic+1]+1  #

                for varobj, v_cdim_type in item["v_info"]:
                    v_name = varobj.name

                    # Copy the subset of data values to the new variable:
                    if v_cdim_type == 10:  # only cdim
                        outp_ds.groups[g_name][v_name][ibc:iec] = (
                                                ds.groups[g_name][v_name][...])
                    elif v_cdim_type%10 == 0:  # cdim first
                        outp_ds.groups[g_name][v_name][ibc:iec,...] = (
                                                ds.groups[g_name][v_name][...])
                    elif v_cdim_type//10 == v_cdim_type%10+1:  # cdim last
                        outp_ds.groups[g_name][v_name][...,ibc:iec] = (
                                                ds.groups[g_name][v_name][...])
                    else:
                        raise ValueError(
                                     f"unimplemented v_cdim_type {v_cdim_type}")

        # Properly re-determine any time coverage global attributes:
        curr_g_atts = outp_ds.ncattrs()
        if "ctime_coverage_start_s" in curr_g_atts:
            outp_ds.setncattr("ctime_coverage_start_s",
                              ds_l[0].ctime_coverage_start_s)
        if "ctime_coverage_end_s" in curr_g_atts:
            outp_ds.setncattr("ctime_coverage_end_s",
                              ds_l[-1].ctime_coverage_end_s)
        if "UTC_coverage_start" in curr_g_atts:
            outp_ds.setncattr("UTC_coverage_start",
                              ds_l[0].UTC_coverage_start)
        if "UTC_coverage_end" in curr_g_atts:
            outp_ds.setncattr("UTC_coverage_end",
                              ds_l[-1].UTC_coverage_end)

    return str(output_fPath)


def dedup_and_concat_along_tdim(raw_fpath_l, output_Path,
                                concatenation_dim_name, dedup_var_info,
                                leap_s_info=None):
    """
    Concatenate multiple NetCDF-format files along a given named time dimension,
     while de-duplicating (with respect to time) the input file stream.

    A single new composite NetCDF-format file (with a name based on the output
     path and the input filepath prefix) results from this operation.

    Input Parameters
    ----------
    raw_fpath_l : list of str
        List of filepaths of the NetCDF-format chunkfiles to be concatenated.
    output_Path : pathlib Path object
        Directory of the output NetCDF-format file.
    concatenation_dim_name : str
        Name of the time dimension (in the input files) to concatenate along
    dedup_var_info : tuple
           (var_name, group_with_var, dedup_v_to_s_factor)
        Name of the (1-D time) variable to use for de-duplication determination,
         the name of the data group that it is in (or None, if no group), and
         the factor to multiply the var values by to get SI seconds (or None, if
         var already has units of SI seconds).
    leap_s_info : dict
        Leap-second reference dictionary
    """

    if len(raw_fpath_l) < 1:
        raise ValueError("no chunk filepaths were provided")

    sorted_fpath_l = sorted(raw_fpath_l)

    # Open all dataset files, read/concatenate raw de-duplication (time) var:
    dedup_vn, dedup_group_w_var, dedup_v_to_s_factor = dedup_var_info
    has_group = (dedup_group_w_var is not None)
    ds_l = []
    for sorted_fpath in sorted_fpath_l:
        ds_l.append(netCDF4.Dataset(sorted_fpath))

        if has_group:
            tmp = ds_l[-1].groups[dedup_group_w_var].variables[dedup_vn][...]
        else:
            tmp = ds_l[-1].variables[dedup_vn][...]

        if len(ds_l) == 1:  # Instantiate array
            raw_dedup_v = tmp.copy()
        else:  # Append to array
            raw_dedup_v = np.concatenate((raw_dedup_v, tmp), axis=0)

    # Sort raw de-duplication (time) variable into ascending time order:
    sorted_inds_full = np.argsort(raw_dedup_v)
    raw_dedup_v_sorted = raw_dedup_v[sorted_inds_full]

    # Detect duplicates, store results as boolean array:
    if dedup_v_to_s_factor is None:
        dedup_threshold_t = (2.e-5, 1.e-5)
    else:
        dedup_threshold_t = (2.e-5/dedup_v_to_s_factor,
                             1.e-5/dedup_v_to_s_factor)
    dcheck = np.append(np.array(dedup_threshold_t[0]),
                       np.diff(raw_dedup_v_sorted))
    dcheck_bool = (dcheck > dedup_threshold_t[1])  # Is this element unique?
    dedup_unique = np.compress(dcheck_bool, raw_dedup_v_sorted, axis=0)
    n_cdim = np.count_nonzero(dcheck_bool)
    print("{} duplicate frames were detected (all: {}, unique: {})".format(
                             len(raw_dedup_v)-n_cdim, len(raw_dedup_v), n_cdim))

    ds0 = ds_l[0]  # Single out the first dataset for copying global atts, etc.

    #== Determine which items should be concatenated and which should only
    #    be copied.  Additionally, determine what dimensions need to be defined
    #    in the output file.

    copyonly_items_gl = []
    concat_items_gl = []
    ds_dim_info = collections.defaultdict(set)

      # For the "root" group:
    g_copy_d = {"g_name": '/'}
    g_copy_d["v_info"] = []
    g_copy_d["a_info"] = ds0.ncattrs()
    g_concat_d = {"g_name": '/'}
    g_concat_d["v_info"] = []
    for v_name in ds0.variables:
        v_dim_names = ds0.variables[v_name].dimensions
        for idim, dim_name in enumerate(v_dim_names):
            ds_dim_info[dim_name].add((ds0.variables[v_name].shape)[idim])

        varobj = ds0.variables[v_name]
        if concatenation_dim_name in v_dim_names:
            v_cdim_type = (len(v_dim_names)*10+
                                     v_dim_names.index(concatenation_dim_name))
            g_concat_d["v_info"].append((varobj, v_cdim_type))
        else:
            g_copy_d["v_info"].append(varobj)
    copyonly_items_gl.append(g_copy_d)
    concat_items_gl.append(g_concat_d)

      # For any other groups:
    for g_name in ds0.groups:
        if len(ds0.groups[g_name].groups) > 0:
            raise RuntimeError("nested groups are not yet supported for "
                               "concatenation")
        g_copy_d = {"g_name": g_name}
        g_copy_d["v_info"] = []
        g_copy_d["a_info"] = ds0.groups[g_name].ncattrs()
        g_concat_d = {"g_name": g_name}
        g_concat_d["v_info"] = []
        for v_name in ds0.groups[g_name].variables:
            v_dim_names = ds0.groups[g_name].variables[v_name].dimensions
            for idim, dim_name in enumerate(v_dim_names):
                ds_dim_info[dim_name].add(
                            (ds0.groups[g_name].variables[v_name].shape)[idim])

            # "v_info" is a list of 'netCDF4.Variable' objects:
            varobj = ds0.groups[g_name].variables[v_name]
            if concatenation_dim_name in v_dim_names:
                v_cdim_type = (len(v_dim_names)*10+
                                     v_dim_names.index(concatenation_dim_name))
                g_concat_d["v_info"].append((varobj, v_cdim_type))
            else:
                g_copy_d["v_info"].append(varobj)
        copyonly_items_gl.append(g_copy_d)
        concat_items_gl.append(g_concat_d)

    #== Determine the output filepath, and create/write it:

    fntoks = os.path.basename(sorted_fpath_l[0]).split('_')
    fn = "{}_ct{}_ct{}.nc".format('_'.join(fntoks[0:4]),
                          int(dedup_unique[0]*1.e2), int(dedup_unique[-1]*1.e2))
    output_fPath = output_Path / fn

    with netCDF4.Dataset(str(output_fPath), 'w', clobber=True) as outp_ds:

        # Create/set all dimensions:
        for dim_name, dim_vset in ds_dim_info.items():
            if dim_name != concatenation_dim_name:
                outp_ds.createDimension(dim_name, list(dim_vset)[0])
        outp_ds.createDimension(concatenation_dim_name, n_cdim)

        # Copy global attributes all at once via dictionary:
        outp_ds.setncatts(ds0.__dict__)

        ## Copy all group attributes and copy-only variables (and their
        ##  attributes):
        for item in copyonly_items_gl:
            g_name = item["g_name"]
            if g_name == '/':  # "root" group variables

                for varobj in item["v_info"]:
                    v_name = varobj.name

                    # Always use shuffle and "zlib" compression:
                    try:
                        outp_ds.createVariable(v_name,
                               varobj.datatype, dimensions=varobj.dimensions,
                               compression="zlib", complevel=4, shuffle=True)
                    except:
                        # Try to use the older 'zlib' argument:
                        outp_ds.createVariable(v_name,
                               varobj.datatype, dimensions=varobj.dimensions,
                               zlib=True, complevel=4, shuffle=True)

                    # Copy variable attributes all at once via dictionary
                    outp_ds[v_name].setncatts(ds0[v_name].__dict__)

                    # Copy the data values to the new variable:
                    outp_ds[v_name][...] = ds0[v_name][...]

            else:  # Other group variables and group attributes
                outp_ds.createGroup(g_name)  # Create each group here

                # Copy any group attributes all at once via dictionary:
                outp_ds.groups[g_name].setncatts(ds0.groups[g_name].__dict__)

                for varobj in item["v_info"]:
                    v_name = varobj.name

                    # Always use shuffle and "zlib" compression:
                    try:
                        outp_ds.groups[g_name].createVariable(v_name,
                               varobj.datatype, dimensions=varobj.dimensions,
                               compression="zlib", complevel=4, shuffle=True)
                    except:
                        # Try to use the older 'zlib' argument:
                        outp_ds.groups[g_name].createVariable(v_name,
                               varobj.datatype, dimensions=varobj.dimensions,
                               zlib=True, complevel=4, shuffle=True)

                    # Copy variable attributes all at once via dictionary
                    outp_ds.groups[g_name][v_name].setncatts(
                                           ds0.groups[g_name][v_name].__dict__)

                    # Copy the data values to the new variable:
                    outp_ds.groups[g_name][v_name][...] = (
                                               ds0.groups[g_name][v_name][...])

        ## Concatenate all relevant variables (copying their attributes):
        for item in concat_items_gl:
            g_name = item["g_name"]
            if g_name == '/':  # "root" group variables

                for varobj, _ in item["v_info"]:
                    v_name = varobj.name
                    v_dimn = varobj.dimensions

                    # Always use shuffle and "zlib" compression:
                    try:
                        outp_ds.createVariable(v_name,
                                 varobj.datatype, dimensions=v_dimn,
                                 compression="zlib", complevel=4, shuffle=True)
                    except:
                        # Try to use the older 'zlib' argument:
                        outp_ds.createVariable(v_name,
                                          varobj.datatype, dimensions=v_dimn,
                                          zlib=True, complevel=4, shuffle=True)

                    # Copy variable attributes all at once via dictionary
                    outp_ds[v_name].setncatts(ds0[v_name].__dict__)

            else:  # Other group variables and group attributes

               for varobj, _ in item["v_info"]:
                    v_name = varobj.name
                    v_dimn = varobj.dimensions

                    # Always use shuffle and "zlib" compression:
                    try:
                        outp_ds.groups[g_name].createVariable(v_name,
                                 varobj.datatype, dimensions=v_dimn,
                                 compression="zlib", complevel=4, shuffle=True)
                    except:
                        # Try to use the older 'zlib' argument:
                        outp_ds.groups[g_name].createVariable(v_name,
                                          varobj.datatype, dimensions=v_dimn,
                                          zlib=True, complevel=4, shuffle=True)

                    # Copy variable attributes all at once via dictionary
                    outp_ds.groups[g_name][v_name].setncatts(
                                           ds0.groups[g_name][v_name].__dict__)

            # Fill the data values of the new output file variables.  Loop
            #  more slowly over the datafiles (than the variables that need
            #  concatenated) in order to avoid unnecessary "drive hopping":
            for varobj, v_cdim_type in item["v_info"]:
                v_name = varobj.name
                raw_v = None

                # Copy the de-duplicated (in time) data values to the new var:
                if v_cdim_type == 10:  # only cdim
                    for ds in ds_l:
                        tmp = ds.groups[g_name].variables[v_name][...]
                        if raw_v is None:  # Instantiate array
                            raw_v = tmp.copy()
                        else:  # Append to array
                            raw_v = np.concatenate((raw_v, tmp), axis=0)
                    raw_v = raw_v[sorted_inds_full]
                    outp_ds.groups[g_name][v_name][...] = np.compress(
                                                     dcheck_bool, raw_v, axis=0)
                elif v_cdim_type%10 == 0:  # cdim first
                    for ds in ds_l:
                        tmp = ds.groups[g_name].variables[v_name][...]
                        if raw_v is None:  # Instantiate array
                            raw_v = tmp.copy()
                        else:  # Append to array
                            raw_v = np.concatenate((raw_v, tmp), axis=0)
                    raw_v = raw_v[sorted_inds_full,...]
                    outp_ds.groups[g_name][v_name][...] = np.compress(
                                                     dcheck_bool, raw_v, axis=0)
                elif v_cdim_type//10 == v_cdim_type%10+1:  # cdim last
                    for ds in ds_l:
                        tmp = ds.groups[g_name].variables[v_name][...]
                        if raw_v is None:  # Instantiate array
                            raw_v = tmp.copy()
                        else:  # Append to array
                            raw_v = np.concatenate((raw_v, tmp),
                                                   axis=len(raw_v.shape)-1)
                    raw_v = raw_v[...,sorted_inds_full]
                    outp_ds.groups[g_name][v_name][...] = np.compress(
                                    dcheck_bool, raw_v, axis=len(raw_v.shape)-1)
                else:
                    raise ValueError(f"unimplemented v_cdim_type {v_cdim_type}")

        # Properly re-determine any time coverage global attributes:
        if "ctime" in dedup_vn:
            curr_g_atts = outp_ds.ncattrs()
            if "ctime_coverage_start_s" in curr_g_atts:
                outp_ds.setncattr("ctime_coverage_start_s",
                                  dedup_unique[0])
            if "ctime_coverage_end_s" in curr_g_atts:
                    outp_ds.setncattr("ctime_coverage_end_s",
                                      dedup_unique[-1])

            if "UTC_coverage" in '\t'.join(curr_g_atts):
                if leap_s_info is not None:
                    UTC_DT, _ = ctime_to_UTC_DT([dedup_unique[0],
                                           dedup_unique[-1]], 's', leap_s_info)
                    UTC_sttrep_t = (UTC_DT[0].strftime("%Y-%m-%dT%H:%M:%S.%f"),
                                    UTC_DT[1].strftime("%Y-%m-%dT%H:%M:%S.%f"))
                else:
                    UTC_sttrep_t = ("unknown", "unknown")

                if "UTC_coverage_start" in curr_g_atts:
                    outp_ds.setncattr("UTC_coverage_start", UTC_sttrep_t[0])
                if "UTC_coverage_end" in curr_g_atts:
                    outp_ds.setncattr("UTC_coverage_end", UTC_sttrep_t[1])

    return str(output_fPath)
