"""
This program requires python version 3.6 or later, and is importable as a
Python module.
"""

  # From the Python standard library:
import os
import sys
import argparse
import datetime

  # From other external Python packages:
import netCDF4
import numpy as np
import shapely.geometry
import pandas as pd
import geopandas
from geocube.api.core import make_geocube
import xarray as xr

  # Custom utilities:
from PREFIRE_PRD_GEN.product_plot_utils import get_ray_angles, \
                                               get_FOV_coords_and_vals, \
                                               read_in_some_granule_data


extcfg_d = {}
extcfg_d["Global"] = {}
extcfg_d["Global"]["output_crs"] = "EPSG:4087"
extcfg_d["Global"]["extent_outpfn_mnk"] = "G"
extcfg_d["Global"]["mpl_ll_extent"] = (
                   [-180., 180., -90., 90.])  # Extent to plot (in polygon_crs)
extcfg_d["Global"]["eff_ll_extent"] = (
      [-180., 180., -90., 90.])  # Extent to extract data from (in polygon_crs)


cfg0_d = {}

cfg0_d["polygon_crs"] = "EPSG:4326"
cfg0_d["plot_only_every_Nth_FOV"] = 1  # Plot all FOVs

# 3413  # EPSG code for "Arctic" polar stereographic
# 4326  # EPSG code for "global" World Geodetic System 1984 (WGS84)
# 4087  # EPSG code ("global" equidistant cylindrical & EPSG:4326/WGS84)
# 3031  # EPSG code for "Antarctic" polar stereographic

cfg0_d["2B-FLX"] = {
        "group": "Flx", "in_var": "olr",
        "specified_channel": None,
        "SRF_fpath": None,
        "out_var": "olr", "field_name": "OLR",
        "fieldv_pxrng": [90., 330.]}  # [W/m^2]
cfg0_d["2B-SFC"] = {
        "group": "Sfc", "in_var": "sfc_spectral_emis",
        "specified_channel": ([13, 21], [11, 20]),
        "SRF_fpath": None,
        "out_var": "sfc_sp_emis", "field_name": "sfc_emis",
        "fieldv_pxrng": ([0.920, 0.995], [0.860, 0.990])}  # [-]
cfg0_d["2B-ATM"] = {
        "group": "Atm", "in_var": "cwv",
        "specified_channel": None,
        "SRF_fpath": None,
        "out_var": "column_wv", "field_name": "column_wv",
        "fieldv_pxrng": [0., 70.]}  # [mm]


#------------------------------------------------------------------------------
def granule_info(p_fpath):
    """Determines top-level info about the given granule."""

    try:
        with netCDF4.Dataset(p_fpath, 'r') as nc_ds:
            sat_num = int(nc_ds.sensor_ID[5])
            granule_IDstr = nc_ds.granule_ID

            granule_type = nc_ds.granule_type
            if granule_type == "1B-NLRAD":
                granule_type = "1B-RAD"
    except:
        msg = f"Input file ({p_fpath}) not found or unreadable"
        raise RuntimeError(msg)

    info_d = {"sat_num": sat_num, "granule_IDstr": granule_IDstr,
              "granule_type": granule_type, "product_item": granule_type}

    return (p_fpath, info_d)


#------------------------------------------------------------------------------
def prepare_polygons_to_plot(ispi, igf, iat_t, product_field_dl, geometry_dl,
                             dims_dl, cfg_d, SRF_fpath=None):
    """Prepare polygons to plot; return a (GeoDataFrame, list, list) tuple."""

    pcfg_d0 = cfg_d[product_field_dl[igf]["product_item"]]
    if pcfg_d0["field_name"] == "brightness_T":
        SRF_d = retrieve_SRF_info(pcfg_d0["SRF_fpath"])
    else:
        SRF_d = None
    
    UTC_beg_l = []
    UTC_end_l = []
    coords_FOV = []
    vals_FOV = []
    dims, geometry, product_field = (dims_dl[igf], geometry_dl[igf],
                                     product_field_dl[igf])

    tmp_coords_FOV, tmp_vals_FOV, UTC_t = get_FOV_coords_and_vals(
                      cfg_d, dims, geometry, product_field, iat_t, ispi, SRF_d)

    coords_FOV.extend(tmp_coords_FOV)
    vals_FOV.extend(tmp_vals_FOV)
    UTC_beg_l.append(UTC_t[0])
    UTC_end_l.append(UTC_t[1])

    df = pd.DataFrame({"geometry": coords_FOV,
                       "plot_field": np.array(vals_FOV, dtype="float32")})
    df["geometry"] = df["geometry"].apply(shapely.geometry.Polygon)
    gdf = geopandas.GeoDataFrame(df, crs=cfg_d["polygon_crs"])

    return (gdf, UTC_beg_l, UTC_end_l)


#------------------------------------------------------------------------------
def plot_product(p_fpath, outp_path, cfg_d):
    """Plot product."""

    # Input parameters:
    #   p_fpath  :  Filepath of product datafile
    #   outp_path  :  path (preferably an absolute one) to write the output
    #                 file(s) to

    # Determine granule info:
    product_fpath_t = granule_info(p_fpath)

    if product_fpath_t[1]["granule_type"] not in cfg_d:
        return  # Nothing more to do here

    # Read in some granule data:
    product_field_dl, geometry_dl, dims_dl = read_in_some_granule_data(
                                                      [product_fpath_t], cfg_d)

    pcfg_d = cfg_d[product_fpath_t[1]["product_item"]]
    pfdl0 = product_field_dl[0]

    try:
        n_spi = len(pfdl0["sp_idx"])
    except:
        n_spi = 1  # scalar sp_idx

      # During the mission we expect <= 6.7 FOV overlaps, so a value of 7 for
      #  this parameter ensures that we do not alias overlapping FOVs:
    iat_stride = 7

    # Bounding box for all polygons:
    gc_bbox = (-180., -90., 180., 90.)
    gc_geom = shapely.geometry.mapping(shapely.geometry.box(*gc_bbox))
    gc_geom["crs"] = {"properties": {"name": cfg_d["polygon_crs"]}}

    gc_fv = -9.985e3  # Wraps to uint8=255
    gc_fv_test = gc_fv+1.

    outp_fpath_l = []
    for ispi in range(n_spi):  # spectral channels

        for igf in range(0, len(dims_dl)):
            for iat_ib in range(0, iat_stride):
                gdf, UTC_beg_l, UTC_end_l = prepare_polygons_to_plot(ispi, igf,
                           (iat_ib, iat_stride), product_field_dl, geometry_dl,
                                                                dims_dl, cfg_d)

                # Rasterize this set of polygons (converting to the output crs):
                gcube = make_geocube(vector_data=gdf, geom=gc_geom,
                                     measurements=["plot_field"], fill=gc_fv,
                                     resolution=(-2.226389815556e3,
                                                 2.226389815556e3),
                                     output_crs=cfg_d["output_crs"])
                
                if iat_ib == 0 and igf == 0:  # Instantiate and initialize
                    sums = np.ma.zeros(gcube.plot_field.values.shape,
                                       dtype="float32")
                    counts = np.zeros(gcube.plot_field.values.shape,
                                      dtype="uint8")

                # Which indices are not masked?  Add their contributions to the
                #  running sum and count:
                gc_inds = (gcube.plot_field > gc_fv_test)
                sums[gc_inds] += gcube.plot_field.values[gc_inds]
                counts[gc_inds] += 1

        sums[counts == 0] = np.ma.masked  # Mask all raster elems with count=0

        # Which raster indices have count > 0?  Divide by the counts to obtain
        #  average raster values for those indices:
        valid_inds = np.nonzero(counts)
        sums[valid_inds] /= counts[valid_inds]
        print ("raw_stats", sums.min(), sums.mean(), sums.max())

        # Rescale raster values to the range 0 to 254, partly determined by
        #  field-dependent min and max values:
        fv_pxrng_min, fv_pxrng_max = tuple(product_field_dl[0]["fv_pxr"])
        sums[sums < fv_pxrng_min] = fv_pxrng_min
        sums[sums > fv_pxrng_max] = (fv_pxrng_max)
        sums = ((sums-fv_pxrng_min)/(fv_pxrng_max-fv_pxrng_min)*254.)

        # Finalize raster values by casting to a datatype of "uint8", and
        #  filling all masked values with 255:
        avg_field_vals = sums.filled(fill_value=255).astype("uint8")

        # Construct an xarray.Dataset using the geocube coordinates:
        ds = xr.Dataset(data_vars={"band1": (('y', 'x'), avg_field_vals)},
                        coords = {'y': gcube.coords['y'].data,
                                  'x': gcube.coords['x'].data})

        # Copy the "band1" DataArray, and add nodata and crs info to it:
        outp_da = ds.band1.rio.set_crs(cfg_d["output_crs"])
        outp_da = outp_da.rio.write_crs(cfg_d["output_crs"], inplace=True)
        outp_da = outp_da.rio.set_nodata(255, inplace=True)
        outp_da = outp_da.rio.write_nodata(255, inplace=True)

        # Construct the output filename and write the raster to file:
        fn = os.path.basename(p_fpath)
        output_fpath = os.path.join(outp_path, "{}.{}{:02d}.tif".format(fn,
                                             cfg_d["extent_outpfn_mnk"], ispi))
        outp_da.rio.to_raster(output_fpath, driver="COG", compress="deflate")
        outp_fpath_l.append(output_fpath)

    return outp_fpath_l


#------------------------------------------------------------------------------
def main(p_fpath, outp_path):
    """Main driver."""

    # Input parameters:
    #   p_fpath  :  Filepath of product datafile
    #   outp_path  :  path (preferably an absolute one) to write the output
    #                 file(s) to

    spec_cfg_d = cfg0_d.copy()
    spec_cfg_d.update(extcfg_d["Global"])
    outp_fpath_l = plot_product(p_fpath, outp_path, spec_cfg_d)

    for fpath in outp_fpath_l:
        print(fpath)

    return outp_fpath_l


#------------------------------------------------------------------------------
if __name__ == "__main__":
    # Process arguments:
    arg_description = ("Create representative-visualization images for NASA "
                       "Worldview.")
    arg_parser = argparse.ArgumentParser(description=arg_description)
    arg_parser.add_argument("product_fpath", help=("Absolute filepath of "
                            "the NetCDF-format product data file to use."))

    args = arg_parser.parse_args()

    output_path = os.path.join(os.getcwd(), "plots")

    # Run driver:
    main(args.product_fpath, output_path)
