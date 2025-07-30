"""
Routines for further describing a NetCDF-format product file (e.g., in a 
separate metadata file).

This program requires python version 3.6 or later, and is importable as a
python module.
"""

  # From the Python standard library:
import os
import json
import hashlib
import datetime

  # From other external Python packages:
import netCDF4
import numpy as np

  # Custom utilities:
import PREFIRE_PRD_GEN.filepaths as filepaths


def produce_DAAC_metadata(output_Path, product_fpath, input_file_type,
                          input_file_info, ext_archival_versionID,
                          aux_input_file_info=None):
    """Produce DAAC per-file/granule metadata (UMM-G format)."""

    # Calculate MD5 hash/checksum for this product file:
    sha256_hash = hashlib.sha256()
    with open(product_fpath, "rb") as f:
        # Read and update hash in chunks of 4K
        for byte_block in iter(lambda: f.read(8192), b""):
            sha256_hash.update(byte_block)
    checksum_d = ({ "Algorithm": "SHA-256",
                    "Value": sha256_hash.hexdigest() })
    
    # Calculate product file size:
    file_size_B = os.path.getsize(product_fpath)  # [B]
    if file_size_B >= 1000000000:
        file_size_rep = file_size_B/1000000000
        file_size_unit = "GB"
    elif file_size_B >= 1000000:
        file_size_rep = file_size_B/1000000
        file_size_unit = "MB"
    if file_size_B >= 1000:
        file_size_rep = file_size_B/1000
        file_size_unit = "KB"
    else:
        file_size_rep = file_size_B
        file_size_unit = "B"

    archival_shortname = "PREFIRE_SAT{}_{}".format(input_file_info["SV_IDstr"],
                                            input_file_type["product_moniker"])
    info_d = {}
    now = datetime.datetime.now(datetime.timezone.utc)
    now_strrep = now.isoformat(sep='T', timespec="milliseconds")[:-6]+'Z'
    if aux_input_file_info is not None:
        assoc_nc_prd_fpath = aux_input_file_info["assoc_nc_prd_fpath"]
          # Open the associated NetCDF-format product file:
        nc_ds = netCDF4.Dataset(assoc_nc_prd_fpath, 'r')
        info_d["archival_versionID"] = nc_ds.archival_versionID
        info_d["UTC_coverage_start"] = nc_ds.UTC_coverage_start
        info_d["UTC_coverage_end"] = nc_ds.UTC_coverage_end
        info_d["UTC_of_file_creation"] = f"{nc_ds.UTC_of_file_creation[:-3]}Z"
        info_d["platform"] = [
                        {"ShortName": f"PREFIRE-SAT{nc_ds.spacecraft_ID[-1]}"}]
        nc_ds.close()  # Close the NetCDF_format product file

        if input_file_type["product_moniker"] == "2B-SFC":
            if ".G00" in product_fpath:
                archival_shortname += "_COG11um"
            elif ".G01" in product_fpath:
                archival_shortname += "_COG19um"
            else:
                raise RuntimeError("unknown coG type")
        else:
            archival_shortname += "_COG"
    elif "NetCDF-format file" in input_file_type:
          # Open the NetCDF-format product file:
        nc_ds = netCDF4.Dataset(product_fpath, 'r')
        info_d["archival_versionID"] = nc_ds.archival_versionID
        try:
            info_d["lat_vertices"] = (nc_ds.groups["Geometry"].
                                             variables["vertex_latitude"][...])
            info_d["lon_vertices"] = (nc_ds.groups["Geometry"].
                                            variables["vertex_longitude"][...])
            info_d["n_atrack"] = nc_ds.dimensions["atrack"].size
            info_d["n_xtrack"] = nc_ds.dimensions["xtrack"].size
        except:
            pass
        info_d["UTC_coverage_start"] = nc_ds.UTC_coverage_start
        info_d["UTC_coverage_end"] = nc_ds.UTC_coverage_end
        info_d["UTC_of_file_creation"] = f"{nc_ds.UTC_of_file_creation[:-3]}Z"
        if input_file_type["product_moniker"] == "AUX-SAT":
            info_d["platform"] = [{
                         "ShortName": f"PREFIRE-SAT{nc_ds.spacecraft_ID[-1]}"},
                                  {"ShortName": "NOAA-20",
                                   "Instruments": [{"ShortName": "VIIRS"}]},
                                  {"ShortName": "Suomi-NPP",
                                   "Instruments": [{"ShortName": "VIIRS"}]},
                                  {"ShortName": "GCOM-W1",
                                   "Instruments": [{"ShortName": "AMSR2"}]},
                                  {"ShortName": "DMSP 5D-3/F18",
                                   "Instruments": [{"ShortName": "SSMIS"}]}]
        else:
            info_d["platform"] = [
                        {"ShortName": f"PREFIRE-SAT{nc_ds.spacecraft_ID[-1]}"}]
        nc_ds.close()  # Close the NetCDF_format product file
    else:
        info_d["archival_versionID"] = ext_archival_versionID
        if info_d["archival_versionID"] is None:
            msg = ("Archival version ID string must be provided for the "
                   "non-NetCDF-format input file {}".format(product_fpath))
            raise RuntimeError(msg)
        info_d["UTC_coverage_start"] = input_file_info["UTC_c_start"]
        info_d["UTC_coverage_end"] = input_file_info["UTC_c_end"]
        info_d["UTC_of_file_creation"] = now_strrep
        info_d["platform"] = [
                    {"ShortName": f"PREFIRE-SAT{input_file_info['SV_IDstr']}"}]

    DAAC_pg_md = {}
    product_fname = os.path.basename(product_fpath)

  # CollectionReference:
    archival_vID = "R{}".format(info_d["archival_versionID"])

    DAAC_pg_md["CollectionReference"] = (
           { "ShortName": archival_shortname,
             "Version": archival_vID })

  # DataGranule:
    DAAC_pg_md["DataGranule"] = {}

    DAAC_pg_md["DataGranule"]["ArchiveAndDistributionInformation"] = ( [
                { "Name": product_fname },
                { "Checksum": checksum_d,
                  "Name": product_fname,
                  "Size": file_size_rep,
                  "SizeUnit": file_size_unit} ] )

    DAAC_pg_md["DataGranule"]["DayNightFlag"] = "Unspecified"
    DAAC_pg_md["DataGranule"]["Identifiers"] = ( [
                { "Identifier": product_fname,
                  "IdentifierType": "ProducerGranuleId" }  ] )

    DAAC_pg_md["DataGranule"]["ProductionDateTime"] = (
                                                info_d["UTC_of_file_creation"])

  # GranuleUR:
    DAAC_pg_md["GranuleUR"] = product_fname

  # MetadataSpecification:
    DAAC_pg_md["MetadataSpecification"] = (
           { "Name": "UMM-G",
             "URL": "https://cdn.earthdata.nasa.gov/umm/granule/v1.6.4",
             "Version": "1.6.4" } )

  # Projects:
    DAAC_pg_md["Projects"] = [ { "ShortName": "PREFIRE" } ]

  # Platforms:
    DAAC_pg_md["Platforms"] = info_d["platform"]

  # ProviderDates:
    DAAC_pg_md["ProviderDates"] = ( [ { "Date": now_strrep,
                                        "Type": "Create" } ] )

  # SpatialExtent:

    try:
          #(atrack, xtrack, vertices):
        lat_vertices = info_d["lat_vertices"]
        lon_vertices = info_d["lon_vertices"]

        n_atrack = info_d["n_atrack"]
        n_xtrack = info_d["n_xtrack"]
        granule_edge_step = max(n_atrack//110, 1)  # Every nth granule edge pt

        granedge_poly_n_pts = (len(range(1, n_atrack, granule_edge_step))+
                               len(range(n_atrack-2, 0, -granule_edge_step))+5)
        granedge_poly_lat = np.empty((granedge_poly_n_pts,))
        granedge_poly_lon = np.empty((granedge_poly_n_pts,))

        i = -1
          # Start of granule to end of granule, along the outer edges of the
          #  first scene's footprints:
        isc = 0
        ia = 0  # (isc,0) for start of granule
        i += 1
        granedge_poly_lat[i], granedge_poly_lon[i] = (lat_vertices[ia,isc,0],
                                                      lon_vertices[ia,isc,0])
        for ia in range(1, n_atrack, granule_edge_step):
            # (isc,3) for mid-granule (arbitrary)
            i += 1
            granedge_poly_lat[i], granedge_poly_lon[i] = (
                                lat_vertices[ia,isc,3], lon_vertices[ia,isc,3])
        if ia != n_atrack-1:
            ia = n_atrack-1  # (isc,3) for end of granule
            i += 1
            granedge_poly_lat[i], granedge_poly_lon[i] = (
                                 lat_vertices[ia,isc,3], lon_vertices[ia,isc,3])

          # End of granule to start of granule, along the outer edges of the
          #  last scene's footprints:
        isc = n_xtrack-1
        ia = n_atrack-1  # (isc,2) for end of granule
        i += 1
        granedge_poly_lat[i], granedge_poly_lon[i] = (lat_vertices[ia,isc,2],
                                                      lon_vertices[ia,isc,2])
        for ia in range(n_atrack-2, 0, -granule_edge_step):
            # (isc,1) for mid-granule (arbitrary)
            i += 1
            granedge_poly_lat[i], granedge_poly_lon[i] = (
                                lat_vertices[ia,isc,1], lon_vertices[ia,isc,1])
        if ia != 0:
            ia = 0  # (isc,1) for start of granule
            i += 1
            granedge_poly_lat[i], granedge_poly_lon[i] = (
                                 lat_vertices[ia,isc,1], lon_vertices[ia,isc,1])

        # Close polygon:
        isc = 0
        ia = 0  # (isc,0) for start of granule
        i += 1
        granedge_poly_lat[i], granedge_poly_lon[i] = (lat_vertices[ia,isc,0],
                                                      lon_vertices[ia,isc,0])
        n_pts = i+1

        granedge_poly_latf = [granedge_poly_lat[i] for i in
                                                         range(n_pts-1, -1, -1)]
        granedge_poly_lonf = [granedge_poly_lon[i] for i in
                                                         range(n_pts-1, -1, -1)]

        pt_l = []
        for i in range(n_pts):
            pt_l.append({ "Latitude": granedge_poly_latf[i],
                          "Longitude": granedge_poly_lonf[i] })

        DAAC_pg_md["SpatialExtent"] = ( {
               "HorizontalSpatialDomain": {
                 "Geometry": {
                   "GPolygons": [
                     {
                       "Boundary": {
                         "Points": pt_l } } ] } } } )
    except:
        if aux_input_file_info is not None:  # Worldview image
            DAAC_pg_md["SpatialExtent"] = ( {
                   "HorizontalSpatialDomain": {
                     "Geometry": {
                       "BoundingRectangles": [
                         { "WestBoundingCoordinate": -180,
                           "NorthBoundingCoordinate": 90,
                           "EastBoundingCoordinate": 180,
                           "SouthBoundingCoordinate": -90 } ] } } } )
        else:
            DAAC_pg_md["SpatialExtent"] = ( {
                   "HorizontalSpatialDomain": {
                     "Geometry": {
                       "BoundingRectangles": [
                         { "WestBoundingCoordinate": -180,
                           "NorthBoundingCoordinate": 84,
                           "EastBoundingCoordinate": 180,
                           "SouthBoundingCoordinate": -84 } ] } } } )

  # TemporalExtent:
    tmp = info_d["UTC_coverage_start"]
    UTC_b = "{}Z".format(tmp[:-3])
    tmp = info_d["UTC_coverage_end"]
    UTC_e = "{}Z".format(tmp[:-3])
    DAAC_pg_md["TemporalExtent"] = ( {
           "RangeDateTime": {
                             "BeginningDateTime": UTC_b,
                             "EndingDateTime": UTC_e } } )

    output_fpath = str(output_Path / os.path.basename(product_fpath))
    with open(output_fpath+".cmr.json", 'w') as out_f:
       json_obj = json.dumps(DAAC_pg_md, indent=2)
       out_f.write(json_obj)
       out_f.write('\n')
