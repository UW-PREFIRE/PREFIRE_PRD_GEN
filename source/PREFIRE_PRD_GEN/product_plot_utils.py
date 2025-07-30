"""
This program requires python version 3.6 or later, and is importable as a
Python module.
"""

  # From the Python standard library:
import datetime

  # From other external Python packages:
import numpy as np
import netCDF4
import shapely.geometry


def get_ray_angles(origin_lond, origin_latd, ray_end_lond_l, ray_end_latd_l):
    """Given coords for an origin point and ray ends, compute the ray angles."""
    # origin_l*d_l :  single longitude or latitude values [deg]
    # ray_end_l*d_l :  (array-like) longitude or latitude values [deg]
    o_lnrad, o_ltrad = tuple(np.deg2rad([origin_lond, origin_latd]))
    re_lnrad = np.deg2rad(ray_end_lond_l)
    re_ltrad = np.deg2rad(ray_end_latd_l)

    dlon = re_lnrad-o_lnrad
    y = np.sin(dlon)*np.cos(re_ltrad)
    x = (np.cos(o_ltrad)*np.sin(re_ltrad)-
         np.sin(o_ltrad)*np.cos(re_ltrad)*np.cos(dlon))

    tmp_bearing = np.rad2deg(np.arctan2(y, x))
    bearing = (tmp_bearing+360.)%360.  # [deg]
#        bearing = 360.-bearing  # count degrees CW - remove to make CCW

    return bearing


def retrieve_SRF_info(SRF_fpath):
    """Extract some useful SRF-related information from the SRF datafile."""
    try:
        with netCDF4.Dataset(SRF_fpath, 'r') as nc_ds:
            NEDR = nc_ds["NEDR"][...]
            rad = nc_ds["rad"][...]
            drad_dT = nc_ds["drad_dT"][...]
            T_grid = nc_ds["T_grid"][...]         
    except:
        raise ValueError("A valid value of 'SRF_fpath' is required.")
    SRF_d = {"NEDR": NEDR, "rad": rad, "drad_dT": drad_dT, "T_grid": T_grid}


def get_FOV_coords_and_vals(cfg_d, dims, geometry, product_field, iat_t, ispi,
                            SRF_d):
    """Prepare FOV polygon info: coordinates of vertices, field value."""
    # SRF_d can be None for a field that does not require SRF-related info.
    pcfg_d = cfg_d[product_field["product_item"]]

    almost_180 = 179.999

    AMline_WhBnd = shapely.geometry.LineString([(-almost_180+360., -85.),
                                                (-almost_180+360., 85.)])
    AMline_EhBnd = shapely.geometry.LineString([(almost_180, -85.),
                                                (almost_180, 85.)])

    # Index ranges to process:
    iat_ib, iat_st = iat_t
    iat_ie = dims["n_atrack"]
    ixt_ib, ixt_ie = (0, dims["n_xtrack"])

    coords_FOV = []
    vals_FOV = []
    
    v_lon = geometry["vertex_lon"][iat_ib:iat_ie:iat_st,ixt_ib:ixt_ie,:]
    v_lat = geometry["vertex_lat"][iat_ib:iat_ie:iat_st,ixt_ib:ixt_ie,:]

    v = pcfg_d["out_var"]
    if pcfg_d["in_var"] in ["spectral_radiance"]:
        if len(product_field[v].shape) == 3:
            sp_rad = product_field[v][iat_ib:iat_ie:iat_st,ixt_ib:ixt_ie,ispi]
            rqf = product_field["rqf"][iat_ib:iat_ie:iat_st,ixt_ib:ixt_ie,ispi]
        else:
            sp_rad = product_field[v][iat_ib:iat_ie:iat_st,ixt_ib:ixt_ie]
            rqf = product_field["rqf"][iat_ib:iat_ie:iat_st,ixt_ib:ixt_ie]
        test1 = ( rqf <= 1 )
        test2 = np.logical_not(np.ma.getmaskarray(sp_rad))
        plot_this = np.logical_and(test1, test2)

        if pcfg_d["field_name"] == "brightness_T":
            scene_IDs = [x+1 for x in range(ixt_ib, ixt_ie)]
            values = btemp_ch(sp_rad, product_field["sp_idx"][ispi]+1,
                              scene_IDs, SRF_d)
        else:
            values = sp_rad
    elif pcfg_d["in_var"] in ["spectral_flux", "sfc_spectral_emis"]:
        if len(product_field[v].shape) == 3:
            sp_rad = product_field[v][iat_ib:iat_ie:iat_st,ixt_ib:ixt_ie,ispi]
        else:
            sp_rad = product_field[v][iat_ib:iat_ie:iat_st,ixt_ib:ixt_ie]
        plot_this = np.logical_not(np.ma.getmaskarray(sp_rad))
        values = sp_rad
    else:
        tmp_field = product_field[v][iat_ib:iat_ie:iat_st,ixt_ib:ixt_ie]
        test1 = np.logical_not(np.ma.getmaskarray(tmp_field))

        test2 = None
        if product_field["in_group"] == "Cld":
            test2 = ( product_field["quality_flag"] <= 1 )

        if test2 is None:
            plot_this = test1
        else:
            plot_this = np.logical_and(test1, test2)

        values = tmp_field

    tUTC = geometry["time_UTC_values"][iat_ib:iat_ie:iat_st,:]
    u = tUTC[0,:]
    UTC_beg = datetime.datetime(u[0], u[1], u[2], u[3], u[4], u[5])
    u = tUTC[-1,:]
    UTC_end = datetime.datetime(u[0], u[1], u[2], u[3], u[4], u[5])

    v_lon_gdiff = np.amax(v_lon, axis=2)-np.amin(v_lon, axis=2)
    split_this_cell = np.where(v_lon_gdiff > 180., True, False)

    glon = geometry["lon"][iat_ib:iat_ie:iat_st,ixt_ib:ixt_ie]
    glat = geometry["lat"][iat_ib:iat_ie:iat_st,ixt_ib:ixt_ie]
    inds = np.nonzero((glon >= cfg_d["eff_ll_extent"][0]) &
                      (glon <= cfg_d["eff_ll_extent"][1]) &
                      (glat >= cfg_d["eff_ll_extent"][2]) &
                      (glat <= cfg_d["eff_ll_extent"][3]))

    for irat, irxt in zip(inds[0], inds[1]):
        if plot_this[irat,irxt]:  # Only plot good detectors
            if (irat%cfg_d["plot_only_every_Nth_FOV"] == 0) and \
                                                    split_this_cell[irat,irxt]:
                # Find vertices with negative (western hemisphere)
                #  longitude values, then construct special FOV polygon
                #  in which longitude can be > 180 deg:
                tmp_v_lat = v_lat[irat,irxt,:]
                tmp_v_lon = v_lon[irat,irxt,:]
                indWh = np.nonzero(tmp_v_lon < 0.)
                otmp_v_lon = v_lon[irat,irxt,:].copy()
                otmp_v_lon[indWh] += 360.
                ocoords_FOV = np.stack((otmp_v_lon, tmp_v_lat), axis=1)
                FOVpoly = shapely.geometry.Polygon(ocoords_FOV)
                    
                # Points for eastern edge of eastern hemisphere subpolygon:
                intersection_pts_EhBnd = list(
                                     FOVpoly.intersection(AMline_EhBnd).coords)

                # Assemble vertex coords for eastern hemisphere subpolygon,
                #  then find an ordering of them that describes a valid
                #  polygon:
                ray_end_lond_l = []  # Initialized here to guard against
                ray_end_latd_l = []  #  case where len(indEht) = 0
                indEht = np.nonzero((tmp_v_lon > 0.) &
                                    (tmp_v_lon < almost_180))
                ray_end_lond_l = [x for x in tmp_v_lon[indEht]]
                ray_end_latd_l = [x for x in tmp_v_lat[indEht]]
                n_ip = len(intersection_pts_EhBnd)
                if n_ip == 2:
                    ray_end_lond_l.append(intersection_pts_EhBnd[0][0])
                    ray_end_lond_l.append(intersection_pts_EhBnd[1][0])
                    ray_end_latd_l.append(intersection_pts_EhBnd[0][1])
                    ray_end_latd_l.append(intersection_pts_EhBnd[1][1])
                elif n_ip == 1:
                    ray_end_lond_l.append(intersection_pts_EhBnd[0][0])
                    ray_end_latd_l.append(intersection_pts_EhBnd[0][1])

                if n_ip > 0:
                    ray_end_lond_a = np.array(ray_end_lond_l)
                    ray_end_latd_a = np.array(ray_end_latd_l)
                   
                    # centroid coords of the convex polygon:
                    cx, cy = (ray_end_lond_a.mean(), ray_end_latd_a.mean())
                    ray_angles = get_ray_angles(cx, cy,
                                                ray_end_lond_a, ray_end_latd_a)
                    sort_inds = np.argsort(ray_angles)

                    # Add eastern hemisphere subpolygon to lists:
                    if len(sort_inds) >= 4:
                        coords_FOV.append(np.stack((ray_end_lond_a[sort_inds],
                                           ray_end_latd_a[sort_inds]), axis=1))
                        vals_FOV.append(values[irat,irxt])
                    
                # Points for western edge of western hemisphere subpolygon:
                tmp_intersection_pts_WhBnd = list(
                                     FOVpoly.intersection(AMline_WhBnd).coords)
                intersection_pts_WhBnd = [(x-360., y) for x, y in
                                          tmp_intersection_pts_WhBnd]

                # Assemble vertex coords for western hemisphere subpolygon,
                #  then find an ordering of them that describes a valid
                #  polygon:
                ray_end_lond_l = []  # Initialized here to guard against
                ray_end_latd_l = []  #  case where len(indWht) = 0
                indWht = np.nonzero((tmp_v_lon < 0.) &
                                    (tmp_v_lon > -almost_180))
                ray_end_lond_l = [x for x in tmp_v_lon[indWht]]
                ray_end_latd_l = [x for x in tmp_v_lat[indWht]]
                n_ip = len(intersection_pts_WhBnd)
                if n_ip == 2:
                    ray_end_lond_l.append(intersection_pts_WhBnd[0][0])
                    ray_end_lond_l.append(intersection_pts_WhBnd[1][0])
                    ray_end_latd_l.append(intersection_pts_WhBnd[0][1])
                    ray_end_latd_l.append(intersection_pts_WhBnd[1][1])
                elif n_ip == 1:
                    ray_end_lond_l.append(intersection_pts_WhBnd[0][0])
                    ray_end_latd_l.append(intersection_pts_WhBnd[0][1])

                if n_ip > 0:
                    ray_end_lond_a = np.array(ray_end_lond_l)
                    ray_end_latd_a = np.array(ray_end_latd_l)
                      # centroid coords of the convex polygon:
                    cx, cy = (ray_end_lond_a.mean(), ray_end_latd_a.mean())
                    ray_angles = get_ray_angles(cx, cy,
                                                ray_end_lond_a, ray_end_latd_a)
                    sort_inds = np.argsort(ray_angles)

                    # Add western hemisphere subpolygon to lists:
                    if len(sort_inds) >= 4:
                        coords_FOV.append(np.stack((ray_end_lond_a[sort_inds],
                                           ray_end_latd_a[sort_inds]), axis=1))
                        vals_FOV.append(values[irat,irxt])
            elif irat%cfg_d["plot_only_every_Nth_FOV"] == 0:
                coords_FOV.append(np.stack((v_lon[irat,irxt,:],
                                            v_lat[irat,irxt,:]), axis=1))
                vals_FOV.append(values[irat,irxt])

    return (coords_FOV, vals_FOV, (UTC_beg, UTC_end))


#------------------------------------------------------------------------------
def get_channel_IDstr_wvl(instrument_str, sp_idx, wvl=None):
    """Return a channel's ID string and possibly its wavelength value(s)."""
    if wvl is None:
        try:
            channel_IDstr = ([f"{instrument_str}" for x
                                               in sp_idx])  # array-like sp_idx
        except:
            channel_IDstr = [f"{instrument_str}"]  # scalar sp_idx
        return channel_IDstr
    else:
        if sp_idx is None:
            channel_IDstr = [f"{instrument_str}"]
            wavelength = None
        else:
            try:
                channel_IDstr = ([f"{x+1:02d}_{instrument_str}" for x
                                               in sp_idx])  # array-like sp_idx
                wavelength = wvl[sp_idx]
            except:
                channel_IDstr = (
                         [f"{sp_idx+1:02d}_{instrument_str}"])  # scalar sp_idx
                wavelength = [wvl[sp_idx]]
        return (channel_IDstr, wavelength)


#------------------------------------------------------------------------------
def read_in_some_granule_data(product_fpath_l, cfg_d):
    """Read in some product granule data."""

    product_field_dl = []
    geometry_dl = []
    dims_dl = []
    for product_fpath, info_d in product_fpath_l:
        product_field = {}
        geometry = {}
        dims = {}
        pcfg_d = cfg_d[info_d["granule_type"]]
        product_field["product_item"] = info_d["product_item"]
        sat_num = info_d["sat_num"]
        fv_pxr = None  # Default value
        if pcfg_d["specified_channel"] is not None:
            sp_idx = pcfg_d["specified_channel"][sat_num-1]
            if "fieldv_pxrng" in pcfg_d:
                fv_pxr = pcfg_d["fieldv_pxrng"][sat_num-1]
        else:
            sp_idx = None
            if "fieldv_pxrng" in pcfg_d:
                fv_pxr = pcfg_d["fieldv_pxrng"]

        with netCDF4.Dataset(product_fpath, 'r') as nc_ds:
            dims["n_xtrack"] = nc_ds.dimensions["xtrack"].size
            dims["n_atrack"] = nc_ds.dimensions["atrack"].size

            geometry["lat"] = (
                           nc_ds.groups["Geometry"].variables["latitude"][...])
            geometry["lon"] = (
                          nc_ds.groups["Geometry"].variables["longitude"][...])

            geometry["vertex_lat"] = (nc_ds.groups["Geometry"].
                                      variables["vertex_latitude"][...])
            geometry["vertex_lon"] = (nc_ds.groups["Geometry"].
                                      variables["vertex_longitude"][...])

            geometry["time_UTC_values"] = (nc_ds.groups["Geometry"].
                                           variables["time_UTC_values"][...])

            geometry["sat_num"] = sat_num
            geometry["instrument"] = f"SAT{sat_num:1d}"
            geometry["granule_ID"] = info_d["granule_IDstr"]
            geometry["full_versionID"] = nc_ds.full_versionID

            product_field["granule_type"] = info_d["granule_type"]
            in_var = pcfg_d["in_var"]
            product_field["in_var"] = in_var
            in_group = pcfg_d["group"]
            product_field["in_group"] = in_group
            out_var = pcfg_d["out_var"]
            product_field["out_var"] = out_var
            product_field["field_name"] = pcfg_d["field_name"]
            if in_var == "spectral_radiance":
                wl = nc_ds.groups[in_group].variables["wavelength"][0,:]
                rqf = (nc_ds.groups[in_group].
                       variables["radiance_quality_flag"][...])

                product_field["channel_ID"], product_field["wavelength"] = (
                          get_channel_IDstr_wvl(geometry["instrument"], sp_idx,
                                                wvl=wl))
                product_field["sp_idx"] = sp_idx
                product_field["fv_pxr"] = fv_pxr
                product_field["rqf"] = rqf[:,:,sp_idx]
                geometry["sat_solar_illumination_flag"] = (
                                             nc_ds.groups["Geometry"].
                                   variables["sat_solar_illumination_flag"][:])
                product_field["observation_bitflags"] = (
                                             nc_ds.groups[in_group].
                                          variables["observation_bitflags"][:])
                product_field[out_var] = (nc_ds.groups[in_group].
                                          variables[in_var][:,:,sp_idx])
            else:
                product_field["sp_idx"] = sp_idx
                product_field["fv_pxr"] = fv_pxr
                try:
                    wl = nc_ds.groups[in_group].variables["wavelength"][0,:]
                    product_field["channel_ID"], product_field["wavelength"] = (
                                  get_channel_IDstr_wvl(geometry["instrument"],
                                                               sp_idx, wvl=wl))
                    product_field[out_var] = (nc_ds.groups[in_group].
                                              variables[in_var][:,:,sp_idx])
                except:
                    product_field["channel_ID"] = get_channel_IDstr_wvl(
                                                geometry["instrument"], sp_idx)
                    product_field[out_var] = (nc_ds.groups[in_group].
                                              variables[in_var][...])
                    if in_group == "Cld":
                        product_field["quality_flag"] = (
                                 nc_ds.groups[in_group].
                                     variables["cld_quality_flag"][:,:])

        product_field_dl.append(product_field)
        geometry_dl.append(geometry)
        dims_dl.append(dims)

    return (product_field_dl, geometry_dl, dims_dl)
