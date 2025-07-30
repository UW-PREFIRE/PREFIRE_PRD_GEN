import netCDF4
import numpy as np


def load_all_vars_of_nc4group(group_name, nc_ds,
                              dim_and_numpy_index_range=None):
    """
    Parameters
    ----------
    group_name : str
        Name of the NetCDF-4 group to load
    nc_ds : netCDF4.Dataset
        A netCDF4.Dataset object referencing the input data file.
    dim_and_numpy_index_range : 3-tuple (str, int, int)
        (optional) A 3-tuple containing the dimension name (in 'nc_ds') to
         subset, and start and stop indices (NumPy indexing convention)

    Returns:
    ----------
    Dictionary containing all variables within the group (possibly subset in
     one of their dimensions).
    """

    gdat = {}
    this_group = nc_ds[group_name]

    if dim_and_numpy_index_range:
        dim_name, ib_np, ie_np = dim_and_numpy_index_range
        for v in this_group.variables:
            try:
                i = this_group.variables[v].dimensions.index(dim_name)
                if i == 0:
                    gdat[v] = this_group.variables[v][ib_np:ie_np,...]
                elif i == 1:
                    gdat[v] = this_group.variables[v][:,ib_np:ie_np,...]
                elif i == 2:
                    gdat[v] = this_group.variables[v][:,:,ib_np:ie_np,...]
                elif i == 3:
                    gdat[v] = this_group.variables[v][:,:,:,ib_np:ie_np,...]
                elif i == 4:
                    gdat[v] = this_group.variables[v][:,:,:,:,ib_np:ie_np,...]
                elif i == 5:
                    gdat[v] = this_group.variables[v][:,:,:,:,:,ib_np:ie_np,...]
            except ValueError:
                gdat[v] = this_group.variables[v][...]
    else:
        for v in this_group.variables:
            gdat[v] = this_group.variables[v][...]

    return gdat


def load_all_atts_of_nc4group(group_name, nc_ds):
    """
    Parameters
    ----------
    group_name : str
        Name of the NetCDF-4 group to load
    nc_ds : netCDF4.Dataset
        A netCDF4.Dataset object referencing the input data file.

    Returns:
    ----------
    Dictionary containing all group attributes of the group.
    """

    gdat = {}
    this_group = nc_ds[group_name]

    for att in this_group.ncattrs():
            gdat[att] = this_group.getncattr(att)

    return gdat


def load_all_nc4_global_atts(nc_ds):
    """
    Parameters
    ----------
    nc_ds : netCDF4.Dataset
        A netCDF4.Dataset object referencing the input data file.

    Returns:
    ----------
    Dictionary containing all global attributes of the file.
    """

    gdat = {}

    for att in nc_ds.ncattrs():
        gdat[att] = nc_ds.getncattr(att)

    return gdat


def get_PREFIRE_Lx_field(nc_ds, group_name, field_name,
                         along_track_np_index_range, as_type=None):
    """
    Extract a given PREFIRE Level-1 or Level-2 field, potentially subsetting it
     in the along-track dimension.

    Parameters
    ----------
    nc_ds : netCDF4.Dataset
        A netCDF4.Dataset object referencing the PREFIRE Level-1 or Level-2
         data file.
    group_name : str or None
        The NetCDF group name that the requested field is a member of (set this
         to None if there is no relevant group)
    field_name : str
        The name of the NetCDF field (variable) that is to be extracted
    along_track_np_index_range : 3-tuple (str, int, int) or None
        A 3-tuple containing the dimension name to subset (in this case,
         "atrack") and start and stop indices (NumPy indexing convention) in
         the given granule to process
        If this is None, do not subset.
    as_type : str or None
        (optional) If this is not None, convert (cast) the extracted field to
                    the specified type (e.g., "float32")

    Returns
    -------
    field_values : numpy.ndarray or numpy.ma.MaskedArray
        NumPy array containing the desired field values (possibly subset, and
         possibly typecast)
    """

    if group_name is None:
        group = '/'
    else:
        group = group_name

    if along_track_np_index_range is None:
        must_subset = False
    else:
        dim_name, ib_np, ie_np = along_track_np_index_range
        must_subset = dim_name in nc_ds.groups[group]. \
                                               variables[field_name].dimensions
    must_typecast = as_type is not None

    if must_subset:
        if must_typecast:
            field_values = nc_ds.groups[group]. \
                         variables[field_name][ib_np:ie_np,...].astype(as_type)
        else:
            field_values = nc_ds.groups[group]. \
                                         variables[field_name][ib_np:ie_np,...]
    else:
        if must_typecast:
            field_values = nc_ds.groups[group]. \
                                     variables[field_name][...].astype(as_type)
        else:
            field_values = nc_ds.groups[group].variables[field_name][...]

    return field_values
