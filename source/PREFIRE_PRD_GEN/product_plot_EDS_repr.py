"""
This program requires python version 3.6 or later, and is importable as a
Python module.
"""

  # From the Python standard library:
import os
import sys
import glob
import argparse
import datetime

  # From other external Python packages:
import matplotlib.figure as mplfig
from matplotlib.colors import ListedColormap
from matplotlib.backends.backend_agg import FigureCanvasAgg
import cartopy
import cartopy.crs as ccrs
import netCDF4
import numpy as np
import shapely.geometry
import pandas as pd
import geopandas

  # Custom utilities:
from PREFIRE_PRD_GEN.product_plot_utils import get_ray_angles, \
                                               get_FOV_coords_and_vals, \
                                               read_in_some_granule_data
from PREFIRE_tools.utils.bitflags import get_bit_from_bitflags_v


#---- Plate Carree map projection --------------------------------------------
cfg_d = {}
cfg_d["plot_only_every_Nth_FOV"] = 1  # Plot all FOVs
cfg_d["proj_args"] = {"central_longitude": 0.}

cfg_d["extent_moniker"] = "Global"
cfg_d["mpl_ll_extent"] = [-180., 180., -90., 90.]
cfg_d["eff_ll_extent"] = [-180., 180., -90., 90.]
cfg_d["coastline_cfg"] = {"scale": "110m", "color": "darkgray"}
cfg_d["ocean_cfg"] = {"scale": "110m", "edgecolor": "darkgray",
                      "facecolor": "none"}
cfg_d["lake_cfg"] = {"scale": "110m", "edgecolor": "darkgray",
                     "facecolor": "none"}

cldmask_cmap = ListedColormap(["#006ba4", "#a2c8ec", "#c85200", "#ababab",
                               "#595959"], name="cldmask")
  # Identical to the first 8 (of 10) values in the color-disabled-friendly
  #  "tableau-colorblind10" color cycle:
auxprod_cmap = ListedColormap(["#006ba4", "#ff800e", "#ababab", "#595959",
                               "#5f9ed1", "#c85200", "#898989", "#a2c8ec"],
                               name="auxprod")

cfg_d["1B-RAD"] = {
        "group": "Radiance", "in_var": "spectral_radiance",
        "specified_channel": (13, 13),
        "SRF_fpath": None,
        "out_var": "sp_radiance", "field_name": "radiance",
        "vminmax_to_use": [["cnd0.", "cnd20."], "neither"],
        "colormap": "viridis",
        "cbar_label": f"Radiance [$W m^{{-2}} sr^{{-1}} \mu m^{{-1}}$]",
        "title_fldname_str": "radiance"}
cfg_d["AUX-MET"] = {
        "group": "Aux-Met", "in_var": "merged_surface_type_prelim",
        "specified_channel": None,
        "SRF_fpath": None,
        "out_var": "mstype_prelim", "field_name": "mstype_prelim",
        "vminmax_to_use": [[0.5, 8.5], "neither"],
        "colormap": auxprod_cmap,
        "cbar_label": ("owat | i$_\mathrm{sea}$ | i$_\mathrm{pse}$ | "
                      "i$_\mathrm{lnd}$ | i$_\mathrm{shf}$ | "
                      "s$_\mathrm{lnd}$ | s$_\mathrm{pln}$ | land"),
        "title_fldname_str": "Merged surface type (GEOS-IT/anc; preliminary)"}
cfg_d["AUX-SAT"] = {
        "group": "Aux-Sat", "in_var": "merged_surface_type_final",
        "specified_channel": None,
        "SRF_fpath": None,
        "out_var": "mstype_final", "field_name": "mstype_final",
        "vminmax_to_use": [[0.5, 8.5], "neither"],
        "colormap": auxprod_cmap,
        "cbar_label": ("owat | i$_\mathrm{sea}$ | i$_\mathrm{pse}$ | "
                      "i$_\mathrm{lnd}$ | i$_\mathrm{shf}$ | "
                      "s$_\mathrm{lnd}$ | s$_\mathrm{pln}$ | land"),
        "title_fldname_str": (
                    "Merged surface type (GEOS-IT/anc/JPSS/AMSR/NISE; final)")}
cfg_d["2B-MSK"] = {
        "group": "Msk", "in_var": "cloud_mask",
        "specified_channel": None,
        "SRF_fpath": None,
        "out_var": "cloud_mask", "field_name": "cloud_mask",
        "vminmax_to_use": [[-0.1, 4.], "neither"],
        "colormap": cldmask_cmap,
        "cbar_label": ("nominal | likely | uncertain | likely | nominal\n"
                "         clear                                cloud         "),
        "title_fldname_str": "retrieved cloud mask"}
cfg_d["2B-FLX"] = {
        "group": "Flx", "in_var": "olr",
        "specified_channel": None,
        "SRF_fpath": None,
        "out_var": "olr", "field_name": "OLR",
        "vminmax_to_use": [["cnd0.", "cnd700."], "neither"],
        "colormap": "viridis",
        "cbar_label": f"Radiative flux [$W m^{{-2}}$]",
        "title_fldname_str": "retrieved outgoing longwave radiative flux"}
cfg_d["2B-SFC"] = {
        "group": "Sfc", "in_var": "sfc_spectral_emis",
        "specified_channel": (25, 25),
        "SRF_fpath": None,
        "out_var": "sfc_sp_emis", "field_name": "sfc_emis",
        "vminmax_to_use": [["cms0.5", "cms1."], "both"],
        "colormap": "cividis",
        "cbar_label": f"Emissivity [-]",
        "title_fldname_str": "retrieved surface (skin) emissivity"}
cfg_d["2B-CLD"] = {
        "group": "Cld", "in_var": "cloud_tau",
        "specified_channel": None,
        "SRF_fpath": None,
        "out_var": "cloud_tau", "field_name": "cloud_tau",
        "vminmax_to_use": [["cnd0.", "cnd10."], "max"],
        "colormap": "viridis",
        "cbar_label": f"Optical depth [-]",
        "title_fldname_str": "retrieved cloud optical depth"}
cfg_d["2B-ATM"] = {
        "group": "Atm", "in_var": "cwv",
        "specified_channel": None,
        "SRF_fpath": None,
        "out_var": "column_wv", "field_name": "column_wv",
        "vminmax_to_use": [["cnd0.", None], "neither"],
        "colormap": "viridis",
        "cbar_label": "Liquid equivalent [mm]",
        "title_fldname_str": "retrieved column water vapor"
        }


#------------------------------------------------------------------------------
def determine_plot_vmn_vmx(pcfg_d, gdf):
    """Determines vmn and vmx for plotting purposes."""

    vmn, vmx = pcfg_d["vminmax_to_use"][0]
    try:
        if vmx[0:3] == "cnd":
            dmax = gdf["plot_field"].max()
            vmx = min(dmax, float(vmx[3:]))
        elif vmx[0:3] == "cms":
            dmean = gdf["plot_field"].mean()
            dstd = gdf["plot_field"].std()
            vmx = min(dmean+dstd, float(vmx[3:]))
    except:
        pass
      # This needs to be after the vmx try block
    try:
        if vmn[0:3] == "cnd":
            dmin = gdf["plot_field"].min()
            vmn = max(dmin, float(vmn[3:]))
        elif vmn[0:3] == "cms":
            dmean = gdf["plot_field"].mean()
            dstd = gdf["plot_field"].std()
            vmn = max(dmean-dstd, float(vmn[3:]))
    except:
        pass

    return (vmn, vmx)


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
def prepare_polygons_to_plot(ispi, product_field_dl, geometry_dl, dims_dl,
                             cfg_d, SRF_fpath=None):
    """Prepare polygons to plot; return a (GeoDataFrame, list, list) tuple."""

    if product_field_dl[0]["field_name"] == "brightness_T":
        SRF_d = retrieve_SRF_info(SRF_fpath)
    else:
        SRF_d = None
    
    UTC_beg_l = []
    UTC_end_l = []
    coords_FOV = []
    vals_FOV = []
    for dims, geometry, product_field in zip(dims_dl, geometry_dl,
                                             product_field_dl):
           # (iat_ib, iat_stride)
        iat_t = (0, 1)

        tmp_coords_FOV, tmp_vals_FOV, UTC_t = get_FOV_coords_and_vals(
                      cfg_d, dims, geometry, product_field, iat_t, ispi, SRF_d)

        coords_FOV.extend(tmp_coords_FOV)
        vals_FOV.extend(tmp_vals_FOV)
        UTC_beg_l.append(UTC_t[0])
        UTC_end_l.append(UTC_t[1])

    df = pd.DataFrame({"geometry": coords_FOV, "plot_field": vals_FOV})
    df["geometry"] = df["geometry"].apply(shapely.geometry.Polygon)
    gdf = geopandas.GeoDataFrame(df)

    gdf = gdf.set_crs(ccrs.PlateCarree())

    return (gdf, UTC_beg_l, UTC_end_l)


#------------------------------------------------------------------------------
def main(p_fpath, outp_path, granule_subset, ancillary_Path=None):
    """Main driver."""

    # Input parameters:
    #   p_fpath  :  Filepath of product datafile
    #   outp_path  :  Path (preferably an absolute one) to write the output
    #                 file(s) to
    #   granule_subset  :  (boolean) Is this a subset granule (and a monitoring
    #                                plot should not be attempted)?
    #   ancillary_Path : (optional; pathlib.Path object) Directory containing
    #                                                    ancillary data files
    #                                                    (cartopy shapefiles)

    # Instruct cartopy to obtain "Natural Earth" shapefile data from a local
    #  source:
    cartopy.config["pre_existing_data_dir"] = str(ancillary_Path / "cartopy")

    # Determine granule info:
    product_fpath_t = granule_info(p_fpath)

    if product_fpath_t[1]["granule_type"] not in cfg_d:
        return  # Nothing more to do here

    # Read in some granule data:
    product_field_dl, geometry_dl, dims_dl = read_in_some_granule_data(
                                                      [product_fpath_t], cfg_d)

    pcfg_d = cfg_d[product_field_dl[0]["granule_type"]]
    produce_mon_plot = ((pcfg_d["field_name"] in ["radiance", "brightness_T"]) &
                        (not granule_subset))

    # Set up a Figure and its canvas -- one set for EDS plots, and possibly
    #  another for monitoring purposes:
    fig_size_t = (11., 7.5)
    fig = []
    canvas = []
    ipEDS = 0
    fig.append(mplfig.Figure(figsize=fig_size_t))
    canvas.append(FigureCanvasAgg(fig[ipEDS]))
    if produce_mon_plot:
        ipmon = 1
        fig.append(mplfig.Figure(figsize=fig_size_t))
        canvas.append(FigureCanvasAgg(fig[ipmon]))

    n_spi = 1

    #-------- plotted on a PlateCarree map:
    mproj_crs = ccrs.PlateCarree(
                     central_longitude=cfg_d["proj_args"]["central_longitude"])

    for ispi in range(n_spi):  # spectral channels
        axo = []
        for ii in range(len(fig)):
            fig[ii].clear()
            axo.append(fig[ii].add_subplot(1, 1, 1, projection=mproj_crs))

            axo[ii].coastlines(resolution=cfg_d["coastline_cfg"]["scale"],
                               color=cfg_d["coastline_cfg"]["color"], zorder=0)

        gdf, UTC_beg_l, UTC_end_l = prepare_polygons_to_plot(ispi,
                                 product_field_dl, geometry_dl, dims_dl, cfg_d)

        pcfg_d = cfg_d[product_field_dl[0]["granule_type"]]

        gdf = gdf.to_crs(mproj_crs)

      #--- EDS plot:

        if pcfg_d["vminmax_to_use"] is None:
            if pcfg_d["cbar_label"] is None:
                gdf.plot(ax=axo[ipEDS], column="plot_field", alpha=0.5,
                         cmap=pcfg_d["colormap"], legend=False)
            else:
                gdf.plot(ax=axo[ipEDS], column="plot_field", alpha=0.5,
                         cmap=pcfg_d["colormap"],
                         legend=True, legend_kwds={"shrink": 0.5, "pad": 0.01,
                                                "label": pcfg_d["cbar_label"]})
        else:
            vmn, vmx = determine_plot_vmn_vmx(pcfg_d, gdf)

            ext_arg = pcfg_d["vminmax_to_use"][1]
            if pcfg_d["cbar_label"] is None:
                gdf.plot(ax=axo[ipEDS], column="plot_field", alpha=0.5,
                         cmap=pcfg_d["colormap"],
                         legend=False, vmin=vmn, vmax=vmx)
            else:
                if '|' in pcfg_d["cbar_label"]:
                    gdf.plot(ax=axo[ipEDS], column="plot_field", alpha=0.5,
                             cmap=pcfg_d["colormap"],
                             legend=True, vmin=vmn, vmax=vmx,
                             legend_kwds={"shrink": 0.5, "pad": 0.01,
                                          "ticks": [vmn, vmx], "format": ' ',
                             "label": pcfg_d["cbar_label"], "extend": ext_arg})
                else:
                    gdf.plot(ax=axo[ipEDS], column="plot_field", alpha=0.5,
                             cmap=pcfg_d["colormap"],
                             legend=True, vmin=vmn, vmax=vmx,
                             legend_kwds={"shrink": 0.5, "pad": 0.01,
                             "label": pcfg_d["cbar_label"], "extend": ext_arg})

      #--- monitoring plot:

        if produce_mon_plot:
            if pcfg_d["vminmax_to_use"] is None:
                gdf.plot(ax=axo[ipmon], column="plot_field", alpha=0.2,
                         cmap="copper", legend=False)
            else:
                vmn, vmx = determine_plot_vmn_vmx(pcfg_d, gdf)

                ext_arg = pcfg_d["vminmax_to_use"][1]
                gdf.plot(ax=axo[ipmon], column="plot_field", alpha=0.2,
                         cmap="copper", legend=False, vmin=vmn, vmax=vmx)

            ref = [(0, "green"), (1, "gold"), (2, "crimson")]
            for ifv, c in ref:
                ind = np.nonzero(product_field_dl[0]["rqf"][:,7] == ifv)
                axo[ipmon].scatter(geometry_dl[0]["lon"][ind,7],
                                   geometry_dl[0]["lat"][ind,7],
                                   transform=mproj_crs, s=0.1, c=c)

            Ttr_after_PObS = get_bit_from_bitflags_v(
                                product_field_dl[0]["observation_bitflags"], 0)
            Tpe_small = get_bit_from_bitflags_v(
                                product_field_dl[0]["observation_bitflags"], 1)
            Tpe_large = get_bit_from_bitflags_v(
                                product_field_dl[0]["observation_bitflags"], 2)
            Tgn_orbit = get_bit_from_bitflags_v(
                                product_field_dl[0]["observation_bitflags"], 3)
            csti_moderate = get_bit_from_bitflags_v(
                                product_field_dl[0]["observation_bitflags"], 4)
            csti_large = get_bit_from_bitflags_v(
                                product_field_dl[0]["observation_bitflags"], 5)

            tmp6 = get_bit_from_bitflags_v(
                                product_field_dl[0]["observation_bitflags"], 6)
            tmp7 = get_bit_from_bitflags_v(
                                product_field_dl[0]["observation_bitflags"], 7)
            bustlm_gapinv = tmp6 | tmp7

            tmp8 = get_bit_from_bitflags_v(
                                product_field_dl[0]["observation_bitflags"], 8)
            tmp9 = get_bit_from_bitflags_v(
                                product_field_dl[0]["observation_bitflags"], 9)
            tmp10 = get_bit_from_bitflags_v(
                               product_field_dl[0]["observation_bitflags"], 10)
            busslew_or_pldwarm = tmp8 | (tmp9 | tmp10)

            ref = [(Ttr_after_PObS, "blue", "T tr after PObS"),
                   (Tpe_small, "red", "T pert, small"),
                   (Tpe_large, "lightcoral", "T pert, large"),
                   (Tgn_orbit, "darkorchid", "g.t.n. orbit T change"),
                   (csti_moderate, "lime", "calseq mod time intv"),
                   (csti_large, "forestgreen", "calseq lrg time intv"),
                   (bustlm_gapinv, "black", "bustlm gap, attdet inv"),
                   (busslew_or_pldwarm, "magenta", "bus slew, pld ewarm")]
            for boolv, c, lb in ref:
                axo[ipmon].scatter(geometry_dl[0]["lon"][boolv,3],
                                   geometry_dl[0]["lat"][boolv,3],
                                   transform=mproj_crs, s=0.1, c=c, label=lb)
            axo[ipmon].legend(fontsize=7., markerscale=3.)

            ind = np.nonzero(geometry_dl[0]["sat_solar_illumination_flag"] == 0)
            axo[ipmon].scatter(geometry_dl[0]["lon"][ind,7],
                               geometry_dl[0]["lat"][ind,7],
                             transform=mproj_crs, s=0.05, c="black", alpha=0.1)

        min_UTC = min(UTC_beg_l)
        fn_dtstr_min = str(min_UTC).replace(' ', 'T').split(':')[0]
        max_UTC = max(UTC_end_l)
        fn_dtstr_max = str(max_UTC).replace(' ', 'T').split(':')[0]

        for ii in range(len(fig)):
            axo[ii].set_extent(cfg_d["mpl_ll_extent"], ccrs.PlateCarree())

            axo[ii].add_feature(cartopy.feature.NaturalEarthFeature(
                        "physical", "ocean", scale=cfg_d["ocean_cfg"]["scale"],
                                          edgecolor="none",
                          facecolor=cfg_d["ocean_cfg"]["facecolor"]), zorder=0)
            axo[ii].add_feature(cartopy.feature.NaturalEarthFeature(
                         "physical", "lakes", scale=cfg_d["lake_cfg"]["scale"],
                                     edgecolor=cfg_d["lake_cfg"]["edgecolor"],
                           facecolor=cfg_d["lake_cfg"]["facecolor"]), zorder=0)
            axo[ii].gridlines(draw_labels=False, xlocs=[],
                              ylocs=[-60., 0., 60.], color="black", alpha=0.3)

            if pcfg_d["field_name"] in ["radiance", "TOA_flux",
                                        "brightness_T", "sfc_emis"]:
                axo[ii].set_title(
                       "PREFIRE_{}_{} {} centered at {:.2f} $\mu$m (ch. {})\n"
                       "(granule {}, {}Z to {}Z)".format(
                           product_field_dl[0]["channel_ID"][ispi][3:],
                           product_field_dl[0]["granule_type"],
                           pcfg_d["title_fldname_str"],
                           product_field_dl[0]["wavelength"][ispi],
                           product_field_dl[0]["channel_ID"][ispi][0:2],
                           geometry_dl[0]["granule_ID"], min_UTC, max_UTC))
            else:
                axo[ii].set_title("PREFIRE_{}_{} {}\n"
                       "(granule {}, {}Z to {}Z)".format(
                           product_field_dl[0]["channel_ID"][ispi],
                           product_field_dl[0]["granule_type"],
                           pcfg_d["title_fldname_str"],
                           geometry_dl[0]["granule_ID"], min_UTC, max_UTC))

        fn = os.path.basename(p_fpath)
        output_fpath = []
        output_fpath.append(os.path.join(outp_path, fn+".jpg"))
        if produce_mon_plot:
            output_fpath.append(os.path.join(outp_path, fn[:-2]+"mon.jpg"))

        for ii in range(len(fig)):
            fig[ii].savefig(output_fpath[ii], bbox_inches="tight", dpi=200)


#------------------------------------------------------------------------------
if __name__ == "__main__":
    # Process arguments:
    arg_description = ("Create representative-visualization images for NASA "
                       "Earthdata search.")
    arg_parser = argparse.ArgumentParser(description=arg_description)
    arg_parser.add_argument("product_fpath", help=("Absolute filepath of "
                            "the NetCDF-format product data file to use."))

    args = arg_parser.parse_args()

    output_path = os.path.join(os.getcwd(), "plots")

    # Run driver:
    main(args.product_fpath, output_path)
