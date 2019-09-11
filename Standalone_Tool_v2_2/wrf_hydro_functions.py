# *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
# ** Copyright UCAR (c) 1992 - 2012
# ** University Corporation for Atmospheric Research(UCAR)
# ** National Center for Atmospheric Research(NCAR)
# ** Research Applications Laboratory(RAL)
# ** P.O.Box 3000, Boulder, Colorado, 80307-3000, USA
# ** 2012/1/5 10:12:13
# *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

import os, sys, zipfile, csv, time
import math, cmath
from zipfile import ZipFile, ZipInfo
from arcpy.sa import *
sys.dont_write_bytecode = True

def makeoutncfile(arcpy, raster, outname, outvar, projdir, loglines):
    outfile = projdir + os.sep + outname
    arcpy.RasterToNetCDF_md(raster, outfile, outvar, "", "x", "y")
    loglines.append('    Process: %s completed without error' %outname)
    arcpy.AddMessage(loglines[-1])
    loglines.append('         Output File: %s' %outfile)
    arcpy.AddMessage(loglines[-1])
    return loglines

def getxy(inraster, projdir, loglines):
    import arcgisscripting
    gp = arcgisscripting.create()
    gp.OverWriteOutput = 1

    # Check out any necessary licenses
    gp.CheckOutExtension("spatial")
    loglines.append('    Starting Process: Converting raster to XMap/YMap')
    gp.ScratchWorkspace = projdir

    OutLyr = 'gplayer'
    gp.MakeRasterLayer_management(inraster, OutLyr)

    # Set environments
    gp.outputCoordinateSystem = OutLyr
    gp.snapRaster= OutLyr
    gp.extent = OutLyr
    gp.cellSize = OutLyr

    # Perform map algebra
    lon = os.path.join(projdir, 'lon')
    lat = os.path.join(projdir, 'lat')
    result1 = gp.SingleOutputMapAlgebra("$$XMap", lon)
    result2 = gp.SingleOutputMapAlgebra("$$YMap", lat)

    # Clean up
    gp.delete_management(OutLyr)
    del OutLyr

    loglines.append('    Conversion of input raster to XMap/YMap completed without error.')
    return lon, lat, loglines

def zipUpFolder(arcpy, zipfile, folder, outZipFile, nclist):
    try:
        zip = zipfile.ZipFile(outZipFile, 'w', zipfile.ZIP_DEFLATED)
        zipws(arcpy, zipfile, str(folder), zip, 'CONTENTS_ONLY', nclist)
        zip.close()
    except RuntimeError:
        pass

def zipws(arcpy, zipfile, path, zip, keep, nclist):
    path = os.path.normpath(path)
    for dirpath, dirnames, filenames in os.walk(path):
        for file in filenames:
            #if file.endswith('.nc'):
            if file in nclist:
                try:
                    if keep:
                        zip.write(os.path.join(dirpath, file), os.path.join(os.sep + os.path.join(dirpath, file)[len(path) + len(os.sep):]))
                except Exception as e:
                    arcpy.AddWarning(get_ID_message(86134) % (file, e[0]))

class ZipCompat(ZipFile):
    def __init__(self, *args, **kwargs):
        ZipFile.__init__(self, *args, **kwargs)

    def extract(self, member, path=None):
        if not isinstance(member, ZipInfo):
            member = self.getinfo(member)
        if path is None:
            path = os.getcwd()
        return self._extract_member(member, path)

    def extractall(self, path=None, members=None, pwd=None):
        if members is None:
            members = self.namelist()
        for zipinfo in members:
            self.extract(zipinfo, path)

    def _extract_member(self, member, targetpath):
        if (targetpath[-1:] in (os.path.sep, os.path.altsep)
            and len(os.path.splitdrive(targetpath)[1]) > 1):
            targetpath = targetpath[:-1]
        if member.filename[0] == '/':
            targetpath = os.path.join(targetpath, member.filename[1:])
        else:
            targetpath = os.path.join(targetpath, member.filename)
        targetpath = os.path.normpath(targetpath)
        upperdirs = os.path.dirname(targetpath)
        if upperdirs and not os.path.exists(upperdirs):
            os.makedirs(upperdirs)
        if member.filename[-1] == '/':
            if not os.path.isdir(targetpath):
                os.mkdir(targetpath)
            return targetpath
        target = file(targetpath, "wb")
        try:
            target.write(self.read(member.filename))
        finally:
            target.close()
        return targetpath

def Examine_Outputs(arcpy, in_zip, out_folder):
    """This tool takes the output zip file from the ProcessGeogrid script and creates a raster
     from each output NetCDF file. The Input should be a .zip file that was created using the
     WRF Hydro pre-processing tools.  The Output Folder parameter should be set to a non-existent
     folder location.  The tool will create the folder which will contain the results."""

    arcpy.AddMessage('Beginning to extract WRF routing grids...')

    # Unzip to a known location (make sure no other nc files live here)
    ZipCompat(in_zip).extractall(out_folder)

    deletefiles = []
    for dirpath, dirnames, filenames in os.walk(out_folder):
        for file in filenames:
            if file.endswith('.nc'):
                # Establish an object for reading the input NetCDF file
                inncfile = os.path.join(out_folder, file)
                ncFP = arcpy.NetCDFFileProperties(inncfile)

                # Loop through global variables in NetCDF file to gather projection information
                ncVariableNames = ncFP.getVariablesByDimension('x')
                variablename = ncVariableNames[-1]
                Dimensions = ncFP.getDimensionsByVariable(variablename)
                Y_Dimension = Dimensions[0]
                X_Dimension = Dimensions[1]
                outRasterLayer = inncfile[:-3]
                arcpy.MakeNetCDFRasterLayer_md(inncfile, variablename, X_Dimension, Y_Dimension, outRasterLayer, "", "", "BY_VALUE")
                arcpy.Raster(outRasterLayer).save()                             # Works!
                #arcpy.Raster(outRasterLayer).save(outRasterLayer + '.tif')
                deletefiles.append(inncfile)
                arcpy.AddMessage('File Created: %s' %outRasterLayer)
                arcpy.Delete_management(outRasterLayer)
                del inncfile, ncFP, ncVariableNames, variablename, Dimensions, Y_Dimension, X_Dimension, outRasterLayer
            if file.endswith('.txt'):
                inASCII = os.path.join(out_folder, file)
                outRasterLayer = os.path.basename(inASCII)[:10] + 'txt'
                outRaster = os.path.join(out_folder, outRasterLayer)
                #outRaster = os.path.join(out_folder, outRasterLayer + '.tif')
                arcpy.ASCIIToRaster_conversion(inASCII, outRaster, 'INTEGER')
                deletefiles.append(inASCII)
                arcpy.AddMessage('File Created: %s' %outRaster)
                del inASCII, outRaster
    del file, dirpath, dirnames, filenames

    # Delete the NC and ASCII input files
    for infile in deletefiles:
        os.remove(infile)
    arcpy.AddMessage('Extraction of WRF routing grids completed.')

def georeference_geogrid_file(arcpy, in_nc, Variable):
    """The first step of the process chain, in which the input NetCDF file gets
    georeferenced and projection files are created"""

    # First step: Import and georeference NetCDF file
    loglines = ['Step 1: NetCDF Conversion initiated... (%s)'%Variable]         # Initiate log list for this process
    arcpy.AddMessage(loglines[-1])

    # Initiate dictionaries of GEOGRID projections and parameters
    projdict = {1: 'Lambert Conformal Conic', 2: 'Polar Stereographic', 3: 'Mercator', 6: 'Cylindrical Equidistant'}
    # See http://www.mmm.ucar.edu/wrf/users/docs/user_guide_V3/users_guide_chap3.htm#_Description_of_the_1

    # Establish an object for reading the input NetCDF file
    ncFP = arcpy.NetCDFFileProperties(in_nc)

    # Loop through global variables in NetCDF file to gather projection information
    ncAttributeNames = ncFP.getAttributeNames("")

    # Find out which projection this GEOGRID file is in
    map_pro = ncFP.getAttributeValue("", 'MAP_PROJ')
    loglines.append('    Map Projection: %s' %projdict[map_pro])
    arcpy.AddMessage(loglines[-1])

    # Collect grid corner XY and DX DY for creating ascii raster later
    if 'corner_lats' in ncAttributeNames:
        corner_lat = ncFP.getAttributeValue("", 'corner_lats')      # Note: The values returned are corner points of the mass grid
    if 'corner_lons' in ncAttributeNames:
        corner_lon = ncFP.getAttributeValue("", 'corner_lons')      # Note: The values returned are corner points of the mass grid
    if 'DX' in ncAttributeNames:
        DX = ncFP.getAttributeValue("", 'DX')
    if 'DY' in ncAttributeNames:
        DY = ncFP.getAttributeValue("", 'DY')

    # Collect necessary information to put together the projection file
    if 'TRUELAT1' in ncAttributeNames:
        standard_parallel_1 = ncFP.getAttributeValue("", 'TRUELAT1')
    if 'TRUELAT2' in ncAttributeNames:
        standard_parallel_2 = ncFP.getAttributeValue("", 'TRUELAT2')
    if 'STAND_LON' in ncAttributeNames:
        central_meridian = ncFP.getAttributeValue("", 'STAND_LON')
    if 'POLE_LAT' in ncAttributeNames:
        pole_latitude = ncFP.getAttributeValue("", 'POLE_LAT')
    if 'POLE_LON' in ncAttributeNames:
        pole_longitude = ncFP.getAttributeValue("", 'POLE_LON')
    if 'CEN_LAT' in ncAttributeNames:
        latitude_of_origin = ncFP.getAttributeValue("", 'CEN_LAT')

    # Process: Make NetCDF Raster Layer
    Dimensions = ncFP.getDimensionsByVariable(Variable)
    Y_Dimension = Dimensions[1]
    X_Dimension = Dimensions[2]
    NC_Raster_Layer = "NC_RASTER"
    arcpy.MakeNetCDFRasterLayer_md(in_nc, Variable, X_Dimension, Y_Dimension, NC_Raster_Layer, "", "", "BY_VALUE")

    # Process: Raster to Numpy Array
    data = arcpy.RasterToNumPyArray(NC_Raster_Layer)
    arcpy.Delete_management(NC_Raster_Layer)
    del NC_Raster_Layer

    # Create Projection file with information from NetCDF global attributes
    sr2 = arcpy.SpatialReference()
    if map_pro == 1:
        # Lambert Conformal Conic
        if 'standard_parallel_2' in locals():
            loglines.append('    Using Standard Parallel 2 in Lambert Conformal Conic map projection.')
            arcpy.AddMessage(loglines[-1])
        else:
            # According to http://webhelp.esri.com/arcgisdesktop/9.2/index.cfm?TopicName=Lambert_Conformal_Conic
            standard_parallel_2 = standard_parallel_1
            latitude_of_origin = standard_parallel_1
        Projection_String = ('PROJCS["Lambert_Conformal_Conic",'
                             'GEOGCS["GCS_Sphere",'
                             'DATUM["D_Sphere",'
                             'SPHEROID["Sphere",6370000.0,0.0]],'
                             'PRIMEM["Greenwich",0.0],'
                             'UNIT["Degree",0.0174532925199433]],'
                             'PROJECTION["Lambert_Conformal_Conic"],'
                             'PARAMETER["false_easting",0.0],'
                             'PARAMETER["false_northing",0.0],'
                             'PARAMETER["central_meridian",' + str(central_meridian) + '],'
                             'PARAMETER["standard_parallel_1",' + str(standard_parallel_1) + '],'
                             'PARAMETER["standard_parallel_2",' + str(standard_parallel_2) + '],'
                             'PARAMETER["latitude_of_origin",' + str(latitude_of_origin) + '],'
                             'UNIT["Meter",1.0]]')

    elif map_pro == 2:
        # Polar Stereographic

        # Set up pole latitude
        phi1 = float(standard_parallel_1)

        ### Back out the central_scale_factor (minimum scale factor?) using formula below using Snyder 1987 p.157 (USGS Paper 1395)
        ##phi = math.copysign(float(pole_latitude), float(latitude_of_origin))    # Get the sign right for the pole using sign of CEN_LAT (latitude_of_origin)
        ##central_scale_factor = (1 + (math.sin(math.radians(phi1))*math.sin(math.radians(phi))) + (math.cos(math.radians(float(phi1)))*math.cos(math.radians(phi))))/2

        # Method where central scale factor is k0, Derivation from C. Rollins 2011, equation 1: http://earth-info.nga.mil/GandG/coordsys/polar_stereographic/Polar_Stereo_phi1_from_k0_memo.pdf
        # Using Rollins 2011 to perform central scale factor calculations. For a sphere, the equation collapses to be much  more compact (e=0, k90=1)
        central_scale_factor = (1 + math.sin(math.radians(abs(phi1))))/2                            # Equation for k0, assumes k90 = 1, e=0. This is a sphere, so no flattening

        loglines.append('      Central Scale Factor: %s' %central_scale_factor)
        arcpy.AddMessage(loglines[-1])
        Projection_String = ('PROJCS["Sphere_Stereographic",'
                             'GEOGCS["GCS_Sphere",'
                             'DATUM["D_Sphere",'
                             'SPHEROID["Sphere",6370000.0,0.0]],'
                             'PRIMEM["Greenwich",0.0],'
                             'UNIT["Degree",0.0174532925199433]],'
                             'PROJECTION["Stereographic"],'
                             'PARAMETER["False_Easting",0.0],'
                             'PARAMETER["False_Northing",0.0],'
                             'PARAMETER["Central_Meridian",' + str(central_meridian) + '],'
                             'PARAMETER["Scale_Factor",' + str(central_scale_factor) + '],'
                             'PARAMETER["Latitude_Of_Origin",' + str(standard_parallel_1) + '],'
                             'UNIT["Meter",1.0]]')

    elif map_pro == 3:
        # Mercator Projection
        Projection_String = ('PROJCS["Sphere_Mercator",'
                             'GEOGCS["GCS_Sphere",'
                             'DATUM["D_Sphere",'
                             'SPHEROID["Sphere",6370000.0,0.0]],'
                             'PRIMEM["Greenwich",0.0],'
                             'UNIT["Degree",0.0174532925199433]],'
                             'PROJECTION["Mercator"],'
                             'PARAMETER["False_Easting",0.0],'
                             'PARAMETER["False_Northing",0.0],'
                             'PARAMETER["Central_Meridian",' + str(central_meridian) + '],'
                             'PARAMETER["Standard_Parallel_1",' + str(standard_parallel_1) + '],'
                             'UNIT["Meter",1.0]]')

    elif map_pro == 6:
        # Cylindrical Equidistant (or Rotated Pole)
        if pole_latitude != float(90) or pole_longitude != float(0):
            # if pole_latitude, pole_longitude, or stand_lon are changed from thier default values, the pole is 'rotated'.
            loglines.append('    Cylindrical Equidistant projection with a rotated pole is not currently supported.')
            arcpy.AddMessage(loglines[-1])
            sys.exit(1)
        else:
            # Check units (linear unit not used in this projection).  GCS?
            Projection_String = ('PROJCS["Sphere_Equidistant_Cylindrical",'
                                 'GEOGCS["GCS_Sphere",'
                                 'DATUM["D_Sphere",'
                                 'SPHEROID["Sphere",6370000.0,0.0]],'
                                 'PRIMEM["Greenwich",0.0],'
                                 'UNIT["Degree",0.0174532925199433]],'
                                 'PROJECTION["Equidistant_Cylindrical"],'
                                 'PARAMETER["False_Easting",0.0],'
                                 'PARAMETER["False_Northing",0.0],'
                                 'PARAMETER["Central_Meridian",' + str(central_meridian) + '],'
                                 'PARAMETER["Standard_Parallel_1",' + str(standard_parallel_1) + '],'
                                 'UNIT["Meter",1.0]]')

    sr2.loadFromString(Projection_String)

    # Create a point geometry object from gathered corner point data
    sr1 = arcpy.SpatialReference(104128)                                        # Using EMEP Sphere (6370000m)
    point = arcpy.Point()
    point.X = corner_lon
    point.Y = corner_lat
    pointGeometry = arcpy.PointGeometry(point, sr1)
    projpoint = pointGeometry.projectAs(sr2)   # Optionally add transformation method:

    # Get projected X and Y of the point geometry and adjust so lower left corner becomes lower left center
    point2 = arcpy.Point((projpoint.firstPoint.X - (DX/2)),(projpoint.firstPoint.Y - (DY/2)))   # Adjust by half a grid cell

    # Process: Numpy Array to Raster
    nc_raster = arcpy.NumPyArrayToRaster(data, point2, DX, DY, -9999)
    del data

    # Process: Define Projection
    arcpy.DefineProjection_management(nc_raster, sr2)
    loglines.append('    Step 1 completed without error.')
    arcpy.AddMessage(loglines[-1])
    return nc_raster, sr2, Projection_String, loglines

def domain_shapefile(arcpy, in_raster, out_shp, sr2):
    """This process creates a shapefile that bounds the GEOGRID file domain. This
    requires the Spatial Analyst extension."""

    if arcpy.CheckExtension("Spatial") == "Available":
        arcpy.CheckOutExtension("Spatial")
        from arcpy.sa import *

    arcpy.AddMessage('Step 2: Build constant raster and convert to polygon...')

    # Set environments
    arcpy.env.overwriteOutput = True
    arcpy.env.outputCoordinateSystem = sr2
    descData = arcpy.Describe(in_raster)
    extent = descData.Extent                                                    # Overkill?
    arcpy.env.snapRaster = in_raster
    cellsize = descData.children[0].meanCellHeight
    arcpy.env.cellSize = in_raster

    # Build constant raster
    constantValue = 1
    outConstRaster = CreateConstantRaster(constantValue, "INTEGER", cellsize, extent)

    # Raster to Polygon conversion
    arcpy.RasterToPolygon_conversion(outConstRaster, out_shp, "NO_SIMPLIFY")    #, "VALUE")

    # Clean up
    arcpy.AddMessage('    Finished building GEOGRID domain boundary shapefile: %s.' %out_shp)

##def create_high_res_topogaphy(arcpy, in_raster, hgt_m_raster, cellsize, sr2, projdir):
##    """The second step creates a high resolution topography raster using a hydrologically-
##    corrected elevation dataset (currently either HydroSHEDS or NHDPlusv2)."""
##
##    if arcpy.CheckExtension("Spatial") == "Available":
##        arcpy.CheckOutExtension("Spatial")
##        from arcpy.sa import *
##
##    # Second part of the process
##    loglines = ['Step 2 initiated...']                                          # Initiate log list for this process
##    arcpy.AddMessage(loglines[-1])
##
##    #Get the extent information from raster object
##    arcpy.MakeRasterLayer_management(hgt_m_raster, 'hgt_m_Layer')
##    descData = arcpy.Describe('hgt_m_Layer')
##    extent = descData.Extent
##    arcpy.env.snapRaster = 'hgt_m_Layer'
##
##    # Test to make sure hgt_m is an integer multiple of supplied output resolution
##    cellsize1 = descData.children[0].meanCellHeight
##    loglines.append('    The GEOGRID File resolution is %sm' %str(cellsize1))
##    arcpy.AddMessage(loglines[-1])
##    cellsize2 = (cellsize1/int(cellsize))
##    loglines.append('    The High-resolution dataset will be %sm' %str(cellsize2))
##    arcpy.AddMessage(loglines[-1])
##
##    # List of coordinates from extent and create a polygon geometry object using an array object
##    coordList = [[extent.XMin, extent.YMin], [extent.XMax, extent.YMin], [extent.XMax, extent.YMax], [extent.XMin, extent.YMax], [extent.XMin, extent.YMin]]
##    boundaryPolygon = arcpy.Polygon(arcpy.Array([arcpy.Point(*coords) for coords in coordList]), sr2)
##
##    # Now project the polygon object to the raster catalog spatial reference
##    sr3 = arcpy.Describe(in_raster).spatialReference
##    transform = arcpy.ListTransformations(sr2, sr3, extent)
##    if len(transform) >= 1:
##        loglines.append('    Tranformation: %s' %transform[0])
##        arcpy.AddMessage(loglines[-1])
##        projpoly = boundaryPolygon.projectAs(sr3, transform[0])                     # Should be: u'NAD_1983_To_WGS_1984_1'
##    else:
##        loglines.append('    Tranformation: None')
##        arcpy.AddMessage(loglines[-1])
##        projpoly = boundaryPolygon.projectAs(sr3)
##    extent2 = projpoly.extent
##
##    # Create raster layer from input raster or mosaic dataset
##    MosaicLayer = "MosaicLayer"
##    arcpy.MakeRasterLayer_management(in_raster, MosaicLayer, "#", extent2)
##    loglines.append('    MakeRasterLayer process completed without error.')
##    arcpy.AddMessage(loglines[-1])
##    loglines.append('    The Coarse Grid has %s rows and %s columns.' %(arcpy.Describe(hgt_m_raster).height, arcpy.Describe(hgt_m_raster).width))
##    arcpy.AddMessage(loglines[-1])
##
##    # Now project the polygon object to the raster catalog spatial reference
##    mosprj = os.path.join(projdir, 'mosaicprj')
##    transform = arcpy.ListTransformations(sr2, sr3, extent)
##    loglines.append('    Projecting input elevation data to WRF coordinate system.')
##    arcpy.AddMessage(loglines[-1])
##    if len(transform) >= 1:
##        arcpy.ProjectRaster_management(MosaicLayer, mosprj, sr2, "NEAREST", cellsize2, transform[0])
##    else:
##        arcpy.ProjectRaster_management(MosaicLayer, mosprj, sr2, "NEAREST", cellsize2, "WGS_1984_(ITRF00)_To_NAD_1983")
##    loglines.append('    Finished projecting input elevation data to WRF coordinate system.')
##    arcpy.AddMessage(loglines[-1])
##    loglines.append('    Fine Grid (before ExtractByMask) has %s rows and %s columns.' %(arcpy.Describe(mosprj).height, arcpy.Describe(mosprj).width))
##    arcpy.AddMessage(loglines[-1])
##
##    # Extract By Mask
##    arcpy.env.cellSize = cellsize2
##    mosprj2 = ExtractByMask(mosprj, hgt_m_raster)                               # Why is this in here? To thin the raster down from the projected raster.
##    arcpy.Delete_management(mosprj)
##    mosprj2.save(os.path.join(projdir, 'mosaicprj'))
##
##    # Check that the number of rows and columns are correct
##    loglines.append('    Fine Grid has %s rows and %s columns.' %(arcpy.Describe(mosprj2).height, arcpy.Describe(mosprj2).width))
##    arcpy.AddMessage(loglines[-1])
##
##    # Clean up
##    arcpy.Delete_management("MosaicLayer")
##    del MosaicLayer
##
##    # Finish
##    loglines.append('    Step 2 completed without error.')
##    arcpy.AddMessage(loglines[-1])
##    return mosprj, cellsize1, cellsize2, loglines

def create_high_res_topogaphy(arcpy, in_raster, hgt_m_raster, cellsize, sr2, projdir):
    """The second step creates a high resolution topography raster using a hydrologically-
    corrected elevation dataset (currently either HydroSHEDS or NHDPlusv2)."""

    if arcpy.CheckExtension("Spatial") == "Available":
        arcpy.CheckOutExtension("Spatial")
        from arcpy.sa import *

    # Second part of the process
    loglines = ['Step 2 initiated...']                                          # Initiate log list for this process
    arcpy.AddMessage(loglines[-1])

    #Get the extent information from raster object
    arcpy.MakeRasterLayer_management(hgt_m_raster, 'hgt_m_Layer')
    descData = arcpy.Describe('hgt_m_Layer')
    extent = descData.Extent
    arcpy.env.snapRaster = 'hgt_m_Layer'                                            # Does this work or does it need to be hgt_m_raster?

    # Test to make sure hgt_m is an integer multiple of supplied output resolution
    cellsize1 = descData.children[0].meanCellHeight
    loglines.append('    The GEOGRID File resolution is %sm' %str(cellsize1))
    arcpy.AddMessage(loglines[-1])
    cellsize2 = (cellsize1/int(cellsize))
    loglines.append('    The High-resolution dataset will be %sm' %str(cellsize2))
    arcpy.AddMessage(loglines[-1])

    # List of coordinates from extent and create a polygon geometry object using an array object
    coordList = [[extent.XMin, extent.YMin], [extent.XMax, extent.YMin], [extent.XMax, extent.YMax], [extent.XMin, extent.YMax], [extent.XMin, extent.YMin]]
    boundaryPolygon = arcpy.Polygon(arcpy.Array([arcpy.Point(*coords) for coords in coordList]), sr2)
    extent1 = arcpy.Extent(coordList[0][0], coordList[0][1], coordList[2][0], coordList[2][1])              # Must do this because extent will change once extent2 gets created (no idea why)

    # Now project the polygon object to the raster catalog spatial reference
    #sr3 = arcpy.Describe(in_raster).spatialReference
    sr3 = sr2
    loglines.append('    Tranformation: faked' )
    transform = arcpy.ListTransformations(sr2, sr3, extent)
    if len(transform) >= 1:
        loglines.append('    Tranformation: %s' %transform[0])
        arcpy.AddMessage(loglines[-1])
        projpoly = boundaryPolygon.projectAs(sr3, transform[0])                     # Should be: u'NAD_1983_To_WGS_1984_1'
    else:
        loglines.append('    Tranformation: None')
        arcpy.AddMessage(loglines[-1])
        projpoly = boundaryPolygon.projectAs(sr3)

    # Create raster layer from input raster or mosaic dataset
    MosaicLayer = "MosaicLayer"
    arcpy.MakeRasterLayer_management(in_raster, MosaicLayer, "#", projpoly.extent)
    loglines.append('    MakeRasterLayer process completed without error.')
    arcpy.AddMessage(loglines[-1])
    loglines.append('    The coarse grid has %s rows and %s columns.' %(arcpy.Describe(hgt_m_raster).height, arcpy.Describe(hgt_m_raster).width))
    arcpy.AddMessage(loglines[-1])
    loglines.append('    The input elevation grid (before reprojection) has %s rows and %s columns.' %(arcpy.Describe(MosaicLayer).height, arcpy.Describe(MosaicLayer).width))
    arcpy.AddMessage(loglines[-1])

    # Set environments to force creation of high-res raster to have exact extent and cellsize needed
    arcpy.env.extent = extent1                                                      # using extent directly doesn't work.
    arcpy.env.outputCoordinateSystem = sr2
    arcpy.env.cellSize = cellsize2
    arcpy.env.snapRaster = hgt_m_raster                                             # Redundant?

    # Now project the polygon object to the raster catalog spatial reference
    mosprj = os.path.join(projdir, 'mosaicprj')
    descData = arcpy.Describe('hgt_m_Layer')
    extent = descData.Extent
    transform = arcpy.ListTransformations(sr2, sr3, extent1)
    loglines.append('    Projecting input elevation data to WRF coordinate system.')
    arcpy.AddMessage(loglines[-1])
    if len(transform) >= 1:
        arcpy.ProjectRaster_management(MosaicLayer, mosprj, sr2, "NEAREST", cellsize2, transform[0])
    else:
        arcpy.ProjectRaster_management(MosaicLayer, mosprj, sr2, "NEAREST", cellsize2, "#","#",sr3)
    loglines.append('      Finished projecting input elevation data to WRF coordinate system.')
    arcpy.AddMessage(loglines[-1])
    loglines.append('    The fine grid (before ExtractByMask) has %s rows and %s columns.' %(arcpy.Describe(mosprj).height, arcpy.Describe(mosprj).width))
    arcpy.AddMessage(loglines[-1])

    # Extract By Mask
    arcpy.env.cellSize = cellsize2
    mosprj2 = ExtractByMask(mosprj, hgt_m_raster)                               # Why is this in here? To thin the raster down from the projected raster.
    arcpy.Delete_management(mosprj)
    mosprj2.save(os.path.join(projdir, 'mosaicprj'))

    # Check that the number of rows and columns are correct
    loglines.append('    Fine Grid has %s rows and %s columns.' %(arcpy.Describe(mosprj2).height, arcpy.Describe(mosprj2).width))
    arcpy.AddMessage(loglines[-1])

    # Clean up
    arcpy.Delete_management("MosaicLayer")
    del MosaicLayer

    # Finish
    loglines.append('    Step 2 completed without error.')
    arcpy.AddMessage(loglines[-1])
    return mosprj, cellsize1, cellsize2, loglines

def create_lat_lon_rasters(arcpy, projdir, mosprj):
    """The third function in the process is to create the latitude and longitude
    rasters that are necessary for running wrf-hydro.  The latitude and longitude
    that are given by the cell values in the resulting output grids are the latitude
    and longitude of the cell center."""

    # Third part of the process
    loglines = ['Step 3 initiated...']
    arcpy.AddMessage(loglines[-1])

    # Execute CreateConstantRaster
    arcpy.env.outputCoordinateSystem = mosprj
    arcpy.env.snapRaster = mosprj
    arcpy.env.extent = mosprj
    arcpy.env.cellSize = mosprj
    OutRas = CreateConstantRaster(1)

    sr1 = arcpy.SpatialReference(104128)         # Using EMEP Sphere (6370000m)
    #sr1 = arcpy.SpatialReference(4326)                                          # Lat and Lon in WGS84
    sr2 = arcpy.Describe(mosprj).spatialReference

    # Project input to output coordinate system
    projraster = os.path.join(projdir, 'mosprj2')
    arcpy.ResetEnvironments()
    arcpy.env.outputCoordinateSystem = arcpy.SpatialReference(104128) # = sr1
    #arcpy.CopyRaster_management(mosprj, projraster)
    arcpy.CopyRaster_management(OutRas, projraster)

    # Create xmap/ymap grids
    lon, lat, loglines = getxy(projraster, projdir, loglines)
    loglines.append('    Latitude and Longitude rasters created.')
    arcpy.AddMessage(loglines[-1])

    arcpy.env.outputCoordinateSystem = mosprj
    arcpy.env.snapRaster = mosprj
    arcpy.env.extent = mosprj
    arcpy.env.cellSize = mosprj

    # Project Raster back to original projection
    CellSize = arcpy.Describe(mosprj).MeanCellHeight
    xout = os.path.join(projdir, 'xoutput')
    yout = os.path.join(projdir, 'youtput')

    # Project the input CRS to the output CRS
    arcpy.ProjectRaster_management(lon, xout, sr2, "NEAREST", CellSize, "#", "#", sr1)
    arcpy.ProjectRaster_management(lat, yout, sr2, "NEAREST", CellSize, "#", "#", sr1)

    ##    transform = arcpy.ListTransformations(sr2, sr1, inraster)
    ##    if len(transform) >= 1:
    ##        arcpy.ProjectRaster_management(lon, xout, sr2, "NEAREST", CellSize, transform[0], "#", sr1)
    ##        arcpy.ProjectRaster_management(lat, yout, sr2, "NEAREST", CellSize, transform[0], "#", sr1)
    ##    else:
    ##        arcpy.ProjectRaster_management(lon, xout, sr2, "NEAREST", CellSize, "WGS_1984_(ITRF00)_To_NAD_1983", "#", sr1)
    ##        arcpy.ProjectRaster_management(lat, yout, sr2, "NEAREST", CellSize, "WGS_1984_(ITRF00)_To_NAD_1983", "#", sr1)

    # Test - extract by  mask
    ##    xout2 = ExtractByMask(xout, mosprj)     # Added 8/2/2014
    ##    yout2 = ExtractByMask(yout, mosprj)     # Added 8/2/2014
    xout2 = ExtractByMask(xout, OutRas)     # Added 8/2/2014
    yout2 = ExtractByMask(yout, OutRas)     # Added 8/2/2014
    arcpy.Delete_management(xout)
    arcpy.Delete_management(yout)

    # Output to NetCDF using Arcpy module
    xoutput = os.path.join(projdir, 'longitude.nc')
    youtput = os.path.join(projdir, 'latitude.nc')
    arcpy.RasterToNetCDF_md(xout2, xoutput, "LONGITUDE", "", "x", "y")
    arcpy.RasterToNetCDF_md(yout2, youtput, "LATITUDE", "", "x", "y")
    arcpy.Delete_management(lon)
    arcpy.Delete_management(lat)
    arcpy.Delete_management(xout2)                  # Added 8/2/2014
    arcpy.Delete_management(yout2)                  # Added 8/2/2014
    arcpy.Delete_management(projraster)

    loglines.append('    Latitude and Longitude NetCDF Files created.')
    arcpy.AddMessage(loglines[-1])
    loglines.append('    Step 3 completed without error.')
    arcpy.AddMessage(loglines[-1])
    return loglines

def Routing_Table(arcpy, projdir, sr2, channelgrid, fdir, Elev, Strahler, loglines):

    """If "Create reach-based routing files?" is selected, this function will create
    the Route_Link.csv table and Streams.shp shapefiles in the output directory."""

    def reorder(ToSeg, Orderlist, Straglers):
        '''This function will navigate through the list of segments until all are accounted
        for. The result is a re-sorted list of which stream segments should be listed
        first. Start with a list of the Order 1 (Strahler) streams. For each link ID,
        the function will find all downstream segments according to the ToSeg dictionary.
        The function will run until there are no more downstream segments. It is important
        to remove circular references or else the function will never complete.'''

        # Fast, order preserving method of uniquifying list (from: http://www.peterbe.com/plog/uniqifiers-benchmark)
        def f5(seq, idfun=None):
           # order preserving
           if idfun is None:
               def idfun(x): return x
           seen = {}
           result = []
           for item in seq:
               marker = idfun(item)
               # in old Python versions: if seen.has_key(marker), but in new ones:
               if marker in seen: continue
               seen[marker] = 1
               result.append(item)
           return result

        # Remove any circular references from the dictionary
        circularlist = [ToSeg.pop(k) for k, v in ToSeg.items() if k==v]             # Ignore self dependencies (Key == Value)

        # Initialize lists
        downlist = list(set(Orderlist))                                             # Elminate duplicates - no need for them. Order of each downlist is not important
        thelist = [1] * len(Orderlist)                                              # Start with a list of ones.
        order = downlist + Straglers                                                # Initiate list with all Order1 streams, then all Order >1 streams but with no contributors (errors)

        # Start navigating downstream. If circular references exist, this will not complete
        count = 1
        while sum(thelist) > 0:
            downlist = [ToSeg.get(x) for x in downlist]                             # Find all downstream segments for the current list order
            order.extend([item for item in downlist if item != None])               # Extend the ordered list by the next group
            thelist = [1 if x > 0 else 0 for x in downlist]                         # Binary list of values (1 == valid value) or None (no downstream segment)
            count += 1

        # Change from ARCID to index
        return f5(order)                                                            # Uses function to preserver order while only keeping unique values

    loglines.append('    Routing table will be created...')
    arcpy.AddMessage(loglines[-1])

    # Get grid information from channelgrid
    descdata = arcpy.Describe(channelgrid)
    extent = descdata.extent
    cellsizeY = descdata.meanCellHeight
    cellsizeX = descdata.meanCellWidth
    del descdata

    # Set output coordinate system environment
    sr1 = arcpy.SpatialReference(104128)                                        # Using EMEP/WRF Sphere (6370000m)
    arcpy.env.outputCoordinateSystem = sr2

    # Output files
    outStreams = os.path.join(projdir, 'Streams.shp')
    RoutingCSV = os.path.join(projdir, 'Route_Link.csv')
    OutFC = os.path.join("in_memory", "Nodes")
    OutFC2 = os.path.join("in_memory", "NodeElev")
    OutFC3 = os.path.join("in_memory", "NodeOrder")
    outRaster = "in_memory/LINKID"

    # Build Stream Features shapefile
    StreamToFeature(channelgrid, fdir, outStreams, "NO_SIMPLIFY")
    loglines.append('        Stream to features step complete.')
    arcpy.AddMessage(loglines[-1])

    # Set environments based on channelgrid
    arcpy.env.snapRaster = channelgrid
    arcpy.env.extent =  channelgrid

    # Create a raster based on the feature IDs that were created in StreamToFeature
    arcpy.FeatureToRaster_conversion(outStreams, 'ARCID', outRaster, channelgrid)                   # Must do this to get "ARCID" field into the raster
    maxValue = arcpy.SearchCursor(outStreams, "", "", "", 'ARCID' + " D").next().getValue('ARCID')  # Gather highest "ARCID" value from field of segment IDs
    maxRasterValue = arcpy.GetRasterProperties_management(outRaster, "MAXIMUM")                     # Gather maximum "ARCID" value from raster
    if maxRasterValue > maxValue:
        whereClause = "VALUE = %s" %maxRasterValue
        outRaster = SetNull(outRaster, outRaster, whereClause)                  # This should eliminate what? NoData values?

    # Create new Feature Class and begin populating it
    arcpy.CreateFeatureclass_management("in_memory", "Nodes", "POINT")
    arcpy.AddField_management(OutFC, "NODE", "LONG")

    # Initiate dictionaries for storing topology information
    From_To = {}                                                                # From_Node/To_Node information
    Nodes = {}                                                                  # Node firstpoint/lastpoint XY information
    NodesLL = {}                                                                # Stores the projected node Lat/Lon information in EMEP Sphere GCS
    Lengths = {}                                                                # Gather the stream feature length
    StrOrder = {}                                                               # Store stream order for each node

    # Enter for loop for each feature/row to gather length and endpoints of stream segments                                                       # Create an array and point object needed to create features
    point = arcpy.Point()
    with arcpy.da.SearchCursor(outStreams, ['SHAPE@', 'ARCID', 'GRID_CODE', 'FROM_NODE', 'TO_NODE']) as rows:                       # Start SearchCursor to look through all linesegments
        for row in rows:
            ID = row[1]                                                         # Get Basin ARCID
            From_To[ID] = [row[3], row[4]]                                      # Store From/To node for each segment
            feat = row[0]                                                       # Create the geometry object 'feat'

            if feat.isMultipart == False:                                       # Make sure that each part is a single part feature
                firstpoint = feat.firstPoint                                    # First point feature geometry
                lastpoint = feat.lastPoint                                      # Last point feature geometry
                Lengths[ID] = feat.length                                       # Store length of the stream segment

                # Gather the X and Y of the top and bottom ends
                for i in firstpoint,lastpoint:                                  # Now put this geometry into new feature class
                    point.X = i.X
                    point.Y = i.Y
                    pointGeometry = arcpy.PointGeometry(point, sr2)
                    projpoint = pointGeometry.projectAs(sr1)                    # Convert to latitude/longitude on the sphere
                    projpoint1 = projpoint.firstPoint
                    if i == firstpoint:                                         # Top Point
                        Nodes[row[3]] = (i.X, i.Y)
                        NodesLL[row[3]] = (projpoint1.X, projpoint1.Y)
                    elif i == lastpoint:                                        # Bottom Point
                        Nodes[row[4]] = (i.X, i.Y)
                        NodesLL[row[4]] = (projpoint1.X, projpoint1.Y)
            else:
                loglines.append('This is a multipart line feature and cannot be handled.')
                arcpy.AddMessage(loglines[-1])
    rows.reset()                                                                # Not sure why this is necessary

    del row, rows, point, ID, feat, firstpoint, lastpoint, pointGeometry, projpoint, projpoint1, i, sr1, sr2
    loglines.append('        Done reading streams layer.')
    arcpy.AddMessage(loglines[-1])

    # Make a point feature class out of the nodes
    IC = arcpy.da.InsertCursor(OutFC, ['SHAPE@', 'NODE'])
    for node in Nodes.keys():                                                   # Now we have to adjust the points that fall outside of the raster edge

        # Adjust X
        if Nodes[node][0] <= extent.XMin:
            Nodes[node] = ((Nodes[node][0] + (cellsizeX/2)), Nodes[node][1])
        elif Nodes[node][0] >= extent.XMax:
            Nodes[node] = ((Nodes[node][0] - (cellsizeX/2)), Nodes[node][1])

        # Adjust Y
        if Nodes[node][1] <= extent.YMin:
            Nodes[node] = (Nodes[node][0], (Nodes[node][1] + (cellsizeY/2)))
        elif Nodes[node][1] >= extent.YMax:
            Nodes[node] = (Nodes[node][0], (Nodes[node][1] - (cellsizeY/2)))

        IC.insertRow([Nodes[node], node])                                       # Insert point and ID information into the point feature class

    del IC, node, Nodes, extent, cellsizeY, cellsizeX
    loglines.append('        Done building Nodes layer with adjustments.')
    arcpy.AddMessage(loglines[-1])

    # Get the elevation values for the nodes feature class
    ExtractValuesToPoints(OutFC, Elev, OutFC2, "NONE", "VALUE_ONLY")
    loglines.append('        Done extracting elevations to points.')
    arcpy.AddMessage(loglines[-1])
    NodeElev = {row[0]: row[1] for row in arcpy.da.SearchCursor(OutFC2, ['NODE', 'RASTERVALU'])}
    loglines.append('        Done reading node elevations.')
    arcpy.AddMessage(loglines[-1])
    arcpy.Delete_management(OutFC2)                                             # Clean up

    # Incorporate Strahler Order
    ExtractValuesToPoints(OutFC, Strahler, OutFC3, "NONE", "VALUE_ONLY")
    with arcpy.da.SearchCursor(OutFC3, ['NODE', 'RASTERVALU']) as rows:
        for row in rows:
            if row[1] <= 0:                                                     # Reclass -9999 values to 1
                order = 1
            else:
                order = row[1]
            StrOrder[row[0]] = order
    loglines.append('        Done reading Strahler stream orders.')
    arcpy.AddMessage(loglines[-1])

    arcpy.Delete_management(OutFC)                                              # Clean up
    arcpy.Delete_management(OutFC3)                                             # Clean up

    # Add stream order into the streams shapefile
    arcpy.AddField_management(outStreams, "Order_", "SHORT")                    # Add field for "Order_"
    with arcpy.da.UpdateCursor(outStreams, ['ARCID', 'Order_']) as rows:        # Start UpdateCursor to add the stream order information
        for row in rows:
            row[1] = StrOrder[From_To[row[0]][0]]
            rows.updateRow(row)

    # Deconstruct from Node space to segment space
    Arc_From = {x:From_To[x][0] for x in From_To}                                       # Build ARCID-keyed topology of The ARCID:FromNode
    From_Arc = {From_To[x][0]:x for x in From_To}                                       # Build Node-keyed topology of The FromNode:ARCID
    From_To2 = {From_To[x][0]:From_To[x][1] for x in From_To}                           # Build Node-keyed topology of The FromNode:ToNode
    Arc_From_To = {item:From_Arc.get(From_To2[Arc_From[item]]) for item in Arc_From}    # Build ARCID-keyed topology of the ARCID:ToARCID
    To_From = {From_To[x][1]:From_To[x][0] for x in From_To}                            # Build Node-keyed topology of The ToNode:FromNode
    Arc_To_From = {item:From_Arc.get(To_From.get(Arc_From[item])) for item in Arc_From} # Build ARCID-keyed topology of the ARCID:FromARCID. ".get()" allows None to be placed in dictionary

    # Get the order of segments according to a simple topological sort
    whereclause = "%s = 1" %arcpy.AddFieldDelimiters(outStreams, "Order_")
    Order1 = [row[0] for row in arcpy.da.SearchCursor(outStreams, 'ARCID', whereclause)] # Generate list of all links of order=1
    Straglers = [row[0] for row in arcpy.da.SearchCursor(outStreams, 'ARCID', whereclause.replace('= 1', '> 1')) if From_To[row[0]][0] is None]  # These are not picked up by the other method
    order = reorder(Arc_From_To, Order1, Straglers)                                     # 'order' variable is a list of LINK IDs that have been reordered according to a simple topological sort

    # Initiate default parameters and build final RouteLink table
    LkHZArea = LkMxH = WeirC = WeirL = OrificC = OrificeA = OrificeE = -9999           # Default parameters
    Qi = 0
    MusK = 3600
    MusX = 0.2
    n = 0.035
    ChSlp = 0.05
    BtmWdth = 5
    type_ = 0
    with open(RoutingCSV, 'wb') as fp:
        a = csv.writer(fp, delimiter=',')
        data = [['link', 'from', 'to', 'start_lon', 'start_lat', 'start_elev', 'type', 'order', 'Qi', 'MusK', 'MusX', 'Length', 'n', 'So', 'ChSlp', 'BtmWdth', 'LkHZArea', 'LkMxH', 'WeirC', 'WeirL', 'OrificC', 'OrificeA', 'OrificeE']]
        for arcid in order:
            fromnode = From_To[arcid][0]                                        # The FROM node from the streams shapefile (used later as a key)
            tonode = From_To[arcid][1]                                          # The TO node from the streams shapefile (used later as a key)
            from_ = Arc_To_From[arcid]                                          # Note that the From in this case is the ARCID of any of the immediately upstream contributing segments
            to = Arc_From_To[arcid]                                             # Note that the To in this case is the ARCID of the immediately downstream segment
            start_lon = round(NodesLL[fromnode][0], 2)                          # Round to 2 digits
            start_lat = round(NodesLL[fromnode][1], 2)                          # Round to 2 digits
            start_elev = round(NodeElev[fromnode], 3)                           # Round to 3 digits
            Length = round(Lengths[arcid], 1)                                   # Round to 1 digit
            drop = NodeElev[fromnode]-NodeElev[tonode]
            # Deal with issue of some segments being assigned higher orders than they should.
            if arcid in Straglers:
                order_ = 1
            else:
                order_ = StrOrder[fromnode]
            if drop < 0:                                                        # Set negative slopes to 0
                drop = 0
            So = round(drop/Length, 3)                                          # Round to 3 digits
            if So == 0.0:
                So == 0.005                                                     # Set minimum slope to be 0.005
            data.append([arcid, from_, to, start_lon, start_lat, start_elev, type_, order_, Qi, MusK, MusX, Length, n, So, ChSlp, BtmWdth, LkHZArea, LkMxH, WeirC, WeirL, OrificC, OrificeA, OrificeE])
        a.writerows(data)
    loglines.append('        Routing Table:\n            %s Lines\n            %s Nodes.' %(len(From_To.keys()), len(NodesLL.keys())))
    arcpy.AddMessage(loglines[-1])
    loglines.append('        Done writing CSV table to disk.')
    arcpy.AddMessage(loglines[-1])

    # Clean up by deletion - Not currently working
    ##    del From_To, NodesLL, NodeElev, StrOrder, Lengths, data, arcid, fp, a
    ##    del Arc_From, From_Arc, From_To2, Arc_From_To, Arc_To_From, To_From
    ##    del from_, to, start_lon, start_lat, start_elev, order_, Length, drop, So, type_, n, ChSlp, Qi, MusK, MusX, BtmWdth, LkHZArea, LkMxH, WeirC, WeirL, OrificC, OrificeA, OrificeE
    ##    del OutFC, OutFC2, OutFC3, RoutingCSV, outStreams, channelgrid, fdir, Elev, Strahler

    loglines.append('    Routing table created without error.')
    arcpy.AddMessage(loglines[-1])
    return outRaster, loglines

def add_reservoirs(arcpy, channelgrid, in_lakes, flac, projdir, fill2, cellsize, sr2, loglines):
    """This function is intended to add reservoirs into the model grid stack, such
    that the channelgrid and lake grids are modified to accomodate reservoirs and
    lakes."""

    # Set lake size threshold
    Threshold = 0.75
    walker = 3

    # Get information about the input domain and set environments
    arcpy.env.cellSize = cellsize
    arcpy.env.extent =  channelgrid
    arcpy.env.outputCoordinateSystem = sr2
    arcpy.env.snapRaster = channelgrid

    # Use extent of channelgrid raster to add a feature layer of lake polygons
    outshp = projdir + os.path.sep + 'in_lakes_clip.shp'
    arcpy.CopyFeatures_management(in_lakes, outshp)

    # Create new area field
    Field1 = 'AREASQKM'
    if Field1 not in [field.name for field in arcpy.ListFields(outshp)]:
        arcpy.AddField_management(outshp, Field1, "FLOAT")
        arcpy.CalculateField_management(outshp, Field1, '!shape.area@squarekilometers!', "PYTHON_9.3")
    field1delim = arcpy.AddFieldDelimiters(outshp, Field1)

    # Use extent of channelgrid raster to add a feature layer of lake polygons
    Field2 = 'FTYPE'
    if Field2 in [field.name for field in arcpy.ListFields(outshp)]:
        Values = ["'LakePond'", "'Reservoir'"]
        field2delim = arcpy.AddFieldDelimiters(outshp, Field2)
        ftypesql = ' AND (%s = %s OR %s = %s)' %(field2delim, Values[0], field2delim, Values[1])
    else:
        ftypesql = ''
    where_clause = """%s >= %s%s""" %(field1delim, Threshold, ftypesql)
    arcpy.MakeFeatureLayer_management(outshp, "Lakeslyr", where_clause)

    # Renumber lakes 1-n (addfield, calculate field)
    newfield = "NEWID"
    arcpy.AddField_management("Lakeslyr", newfield, "LONG")
    expression = 'autoIncrement()'
    code_block = """rec = 0\ndef autoIncrement():\n    global rec\n    pStart = 1\n    pInterval = 1\n    if (rec == 0):\n        rec = pStart\n    else:\n        rec = rec + pInterval\n    return rec"""
    arcpy.CalculateField_management("Lakeslyr", newfield, expression, "PYTHON_9.3", code_block)

    # Create a raster from the lake polygons that matches the channelgrid layer
    outRastername = os.path.join(projdir, "Lakesras")
    outfeatures = os.path.join(projdir, 'outfeatures.shp')
    arcpy.CopyFeatures_management("Lakeslyr", outfeatures)
    arcpy.PolygonToRaster_conversion(outfeatures, newfield, outRastername, "MAXIMUM_AREA")      # This tool requires ArcGIS for Desktop Advanced OR Spatial Analyst Extension
    #arcpy.Featurearcpy.FeatureToRaster_conversion(outfeatures, newfield, outRastername)        # This tool requires only ArcGIS for Desktop Basic, but does not allow a priority field

    # Hack to convert Lakesras to 16bit integer
    outRaster1 = (channelgrid * 0) + Raster(outRastername)                      # Create 16bit ratser for Lakesras out of channelgrid
    outRaster = Con(IsNull(outRaster1)==1, -9999, outRaster1)                   # Convert Null or NoData to -9999

    # Get COMID from lakes feature layer
    COMIDs = {row[0]: row[1] for row in arcpy.da.SearchCursor("Lakeslyr", [newfield, 'COMID'])}     # Currently, COMID must be in the lakes layer

    # Con statement for lakes where channelgrid = -9999
    zonstat = ZonalStatistics(outRastername, "Value", flac, "MAXIMUM")          # Get maximum flow accumulation value for each lake
    lakeacc = Con(outRastername, flac)                                          # Flow accumulation over lakes only
    TestCon = Con(lakeacc == zonstat, outRastername, -9999)                     # Bottom of lake channelgrid pixel = lake number
    NewChannelgrid = Con(IsNull(TestCon) == 1, channelgrid, TestCon)            # This is the new lake-routed channelgrid layer

    # Now march down 3 pixels to get minimum lake elevation
    tolerance = int(float(arcpy.GetRasterProperties_management(channelgrid, 'CELLSIZEX').getOutput(0)) * walker)
    SnapPour = SnapPourPoint(SetNull(TestCon, TestCon, "VALUE = -9999"), flac, tolerance)                          # Snap lake outlet pixel to FlowAccumulation with tolerance
    outTable2 = "in_memory/zonalstat2"
    Sample(fill2, SnapPour, outTable2, "NEAREST")
    min_elevs = {row[1]: row[-1] for row in arcpy.da.SearchCursor(outTable2, "*")}

    # Convert lakes raster to polygons and gather size & elevations
    loglines.append('    Gathering lake parameter information.')
    arcpy.AddMessage(loglines[-1])
    outTable = "in_memory/zonalstat"
    zontable = ZonalStatisticsAsTable(outRastername, 'VALUE', fill2, outTable, "DATA", "MIN_MAX_MEAN")
    areas = {row[0]: row[1] for row in arcpy.da.SearchCursor(outTable, ['VALUE', 'AREA'])}                      # Searchcursor on zonal stats table
    max_elevs = {row[0]: row[1] for row in arcpy.da.SearchCursor(outTable, ['VALUE', 'MAX'])}                   # Searchcursor on zonal stats table

    OrificEs = {x:(min_elevs[x] + ((max_elevs[x] - min_elevs[x])/3)) for x in min_elevs.keys()}             # OrificElevation is 1/3 between the low elevation and max lake elevation
    Elevation = {x:(min_elevs[x] + (((max_elevs[x] - min_elevs[x])/3) * 2)) for x in min_elevs.keys()}      # Elevation is 2/3 between the low elevation and max lake elevation

    #  Gather centroid lat/lons
    out_lake_raster = os.path.join(projdir, "out_lake_raster.shp")
    out_lake_raster_dis = os.path.join(projdir, "out_lake_raster_dissolve.shp")
    arcpy.RasterToPolygon_conversion(outRastername, out_lake_raster, "NO_SIMPLIFY", "VALUE")
    arcpy.Dissolve_management(out_lake_raster, out_lake_raster_dis, "GRIDCODE", "", "MULTI_PART")               # Dissolve to eliminate multipart features

    # Create a point geometry object from gathered lake centroid points
    loglines.append('    Starting to gather lake centroid information.')
    arcpy.AddMessage(loglines[-1])
    sr1 = arcpy.SpatialReference(4326)                                          # GCS_WGS_1984
    point = arcpy.Point()
    cen_lats = {}
    cen_lons = {}
    shapes = {row[0]: row[1] for row in arcpy.da.SearchCursor(out_lake_raster_dis, ['GRIDCODE', 'SHAPE@XY'])}
    for shape in shapes:
        point.X = shapes[shape][0]
        point.Y = shapes[shape][1]
        pointGeometry = arcpy.PointGeometry(point, sr2)
        projpoint = pointGeometry.projectAs(sr1)                                # Optionally add transformation method:
        cen_lats[shape] = projpoint.firstPoint.Y
        cen_lons[shape] = projpoint.firstPoint.X
    loglines.append('    Done gathering lake centroid information.')
    arcpy.AddMessage(loglines[-1])

    # Create Lake parameter file
    LakeCSV = os.path.join(projdir, 'LAKEPARM.TBL')
    OrificC = 0.1
    OrificA = 1.0
    WeirC = 0.4
    WeirL = 0.0
    loglines.append('    Starting to create lake parameter table.')
    arcpy.AddMessage(loglines[-1])
    with open(LakeCSV, 'wb') as fp:
        a = csv.writer(fp, dialect='excel-tab', quoting=csv.QUOTE_NONE)
        a.writerow(['lake', 'LkArea', 'LkMxH', 'WeirC', 'WeirL', 'OrificC', 'OrificeA', 'OrificeE', 'lat', 'long', 'elevation'])
        for lkid in min_elevs.keys():
            lkarea = float(areas[lkid])/float(1000000)                          # Divide by 1M for kilometers^2
            lkmaxelev = max_elevs[lkid]
            baseelev = Elevation[lkid]                                          # Base Elevation is 2/3 betwen 'min' and max lake elevation.
            OrificeE = OrificEs[lkid]                                           # Orifice Elevation is 1/3 between 'min' and max lake elevation.
            #COMID = COMIDs[lkid]
            cen_lat = cen_lats[lkid]
            cen_lon = cen_lons[lkid]
            a.writerow([lkid, lkarea, lkmaxelev, WeirC, WeirL, OrificC, OrificA, OrificeE, cen_lat, cen_lon, baseelev])   #COMID?
    loglines.append('        Lakes Table: %s Lakes' %len(areas.keys()))
    arcpy.AddMessage(loglines[-1])
    loglines.append('        Done writing LAKEPARM.TBL table to disk.')
    arcpy.AddMessage(loglines[-1])

    # Clean up by deletion
    loglines.append('    Lake parameter table created without error.')
    arcpy.AddMessage(loglines[-1])

    # Clean up
    del sr1, point, outRaster1
    arcpy.Delete_management("Lakeslyr")
    arcpy.Delete_management(outshp)
    arcpy.Delete_management(outfeatures)
    arcpy.Delete_management('in_memory')

    outRaster.save(outRastername)
    return arcpy, NewChannelgrid, outRastername, loglines

def sa_functions(arcpy, basin_mask, mosprj, ovroughrtfac_val, retdeprtfac_val, projdir, in_csv, out_zip, threshold, inunits, LU_INDEX, cellsize1, cellsize2, routing, in_lakes):
    """The last major function in the processing chain is to perform the spatial
    analyst functions to hydrologically process the input raster datasets."""

    if arcpy.CheckExtension("Spatial") == "Available":
        arcpy.CheckOutExtension("Spatial")
        from arcpy.sa import *

    # Fourth part of the process
    loglines = ['Step 4 initiated...']
    arcpy.AddMessage(loglines[-1])

    # Set Basin mask attribute to boolean from ArcGIS text
    if basin_mask == 'true':
        bsn_msk = True
        loglines.append('    Channelgrid will be masked to basins.')
        arcpy.AddMessage(loglines[-1])
    elif basin_mask == 'false':
        bsn_msk = False
        loglines.append('    Channelgrid will not be masked to basins.')
        arcpy.AddMessage(loglines[-1])

    # Set routing attribute to boolean from ArcGIS text
    if routing == 'true':
        routing_table = True
        loglines.append('    Reach-based routing files will be created.')
        arcpy.AddMessage(loglines[-1])
    elif routing == 'false':
        routing_table = False
        loglines.append('    Reach-based routing files will not be created.')
        arcpy.AddMessage(loglines[-1])

    walker = 3                                                                  # Number of cells to walk downstream before watershed delineation

    # Check out licenses and set environment variables
    arcpy.CheckOutExtension('spatial')
    arcpy.env.overwriteOutput = True

    # Set environments
    arcpy.MakeRasterLayer_management(mosprj, 'mosaicprj')
    sr2 = arcpy.Describe('mosaicprj').spatialReference
    arcpy.env.outputCoordinateSystem = sr2

    fill = Fill(mosprj)
    #fill = mosprj
    if inunits == 'm':
        fill1 = Float(fill)
        del fill
    elif inunits == 'cm':                                                       # trap for fixing cm to m conversion
        fill1 = fill/float(100)
        del fill

    # Process: Flow Direction
    fdir = FlowDirection(fill1)
    loglines = makeoutncfile(arcpy, fdir, 'flowdirection.nc', 'FLOWDIRECTION', projdir, loglines)

    # Process: Flow Accumulation (intermediate
    flac = FlowAccumulation(fdir, '#', 'FLOAT')
    loglines = makeoutncfile(arcpy, flac, 'flowacc.nc', 'FLOWACC', projdir, loglines)

    # Set NoData Value for topography
    fill2 = Con(IsNull(fill1) == 1, -9999, fill1)
    del fill1
    nodataval = flac.noDataValue
    arcpy.SetRasterProperties_management(fill2, "ELEVATION", "#", "#", [[1, nodataval]])     # Must set NoData away from -Infinityf
    outfile = projdir + os.path.sep + 'topography.nc'
    outfile_org = projdir + os.path.sep + 'topography_org.nc'
    arcpy.RasterToNetCDF_md(fill2, outfile, 'TOPOGRAPHY', "", "x", "y")
    arcpy.RasterToNetCDF_md(mosprj, outfile_org, 'TOPOGRAPHY', "", "x", "y")
    loglines.append('    Process topography.nc completed without error')
    arcpy.AddMessage(loglines[-1])
    loglines.append('         Output File: %s' %outfile)
    loglines.append('         Output File: %s' %outfile_org)
    arcpy.AddMessage(loglines[-1])

    # Create stream channel raster according to threshold
    strm = SetNull(flac, '1', 'VALUE < %s' % threshold)
    channelgrid = Con(IsNull(strm) == 0, 0, -9999)

    # Uncomment these if you want to hardcode the script to use a certain channelgrid raster
    #strm = arcpy.Raster(r'E:\Projects\domains\URG\URGRB_2014_10_30\Experiments\newstrm1')
    #channelgrid = arcpy.Raster(r'E:\Projects\domains\URG\URGRB_2014_10_30\Experiments\newchnlgrd3')

    # Create initial constant raster of -9999
    constraster = CreateConstantRaster(-9999, "INTEGER")

    # Create initial constant raster of value retdeprtfac_val
    retdeprtfac_value = float(retdeprtfac_val)
    inraster2 = CreateConstantRaster(retdeprtfac_value, "FLOAT")
    loglines = makeoutncfile(arcpy, inraster2, 'retdeprtfac.nc', 'RETDEPRTFAC', projdir, loglines)
    del retdeprtfac_value, inraster2

    # Create initial constant raster of ovroughrtfac_val
    ovroughrtfac_value = float(ovroughrtfac_val)
    inraster3 = CreateConstantRaster(ovroughrtfac_value, "FLOAT")
    loglines = makeoutncfile(arcpy, inraster3, 'ovroughrtfac.nc', 'OVROUGHRTFAC', projdir, loglines)
    del ovroughrtfac_value, inraster3

    # Process: Stream Order
    order = StreamOrder(strm, fdir)         # Default = "STRAHLER"
    order2 = Con(IsNull(order) == 1, -9999, order)
    loglines = makeoutncfile(arcpy, order2, 'str_order.nc', 'STREAMORDER', projdir, loglines)

    if routing_table == True:
        linkid, loglines = Routing_Table(arcpy, projdir, sr2, strm, fdir, fill2, order2, loglines)
        loglines = makeoutncfile(arcpy, linkid, 'LINKID.nc', 'LINKID', projdir, loglines)
    del order, order2

    # Find out if forecast points are chosen, then set mask for them
    if in_csv:
        # Make feature layer from CSV
        loglines.append('    Forecast points provided and basins being delineated.')
        arcpy.AddMessage(loglines[-1])
        frxst_layer = 'frxst_layer'
        sr1 = arcpy.SpatialReference(4326)              # GCS_WGS_1984
        #loglines.append('    Spatial Reference: %s' %sr1.exportToString())
        #arcpy.AddMessage(loglines[-1])
        arcpy.MakeXYEventLayer_management(in_csv, 'LON', 'LAT', frxst_layer, sr1)
        tolerance = int(float(arcpy.GetRasterProperties_management('mosaicprj', 'CELLSIZEX').getOutput(0)) * walker)
        tolerance1 = int(float(arcpy.GetRasterProperties_management('mosaicprj', 'CELLSIZEX').getOutput(0)))
        frxst_raster = SnapPourPoint(frxst_layer, flac, tolerance1, 'FID')
        frxst_raster2 = Con(IsNull(frxst_raster) == 0, frxst_raster, -9999)
        loglines = makeoutncfile(arcpy, frxst_raster2, 'frxst_pts.nc', 'frxst_pts', projdir, loglines)
        SnapPour = SnapPourPoint(frxst_layer, flac, tolerance)
        arcpy.Delete_management(frxst_layer)
        arcpy.Delete_management(frxst_raster2)

        # Delineate above points
        outWatershed = Watershed(fdir, SnapPour, 'VALUE')
        watershedraster = os.path.join(projdir, 'watersheds')
        outWatershed.save(watershedraster)
        outWatershed2 = Con(IsNull(outWatershed) == 0, outWatershed, -9999)
        loglines = makeoutncfile(arcpy, outWatershed2, 'frxst_basns.nc', 'basn_msk', projdir, loglines)

        # Groundwater basins are for now the same as the forecast basins
        loglines = makeoutncfile(arcpy, outWatershed2, 'gw_basns.nc', 'basn_msk', projdir, loglines)

        # Process: Resample gw_basns grid to a lower resolution
        gw_basns_geogrid = os.path.join(projdir, 'watersheds2')
        outASCII = os.path.join(projdir, 'gw_basns_geogrid.txt')
        arcpy.env.cellSize = cellsize1                                          # Set cellsize environment to coarse grid
        arcpy.Resample_management(outWatershed2, gw_basns_geogrid, cellsize1, "MAJORITY")
        arcpy.RasterToASCII_conversion(gw_basns_geogrid, outASCII)
        arcpy.env.cellSize = cellsize2                                          # Set cellsize environment back to fine grid
        loglines.append('    Process: gw_basns_geogrid.txt completed without error')
        arcpy.AddMessage(loglines[-1])
        loglines.append('         Output File: %s' %outASCII)
        arcpy.AddMessage(loglines[-1])

        # Set mask for future raster output
        if bsn_msk == True:
            channelgrid2 = Con(outWatershed2 >= 0, IsNull(strm), Con(IsNull(strm), 1, -1))  # Converts channelgrid values inside basins to 0, outside to -1
            channelgrid = Con(channelgrid2 == 1, -9999, channelgrid2)

    else:
        # Handle the change if input forecast points are not provided
        loglines = makeoutncfile(arcpy, constraster, 'frxst_pts.nc', 'frxst_pts', projdir, loglines)
        loglines.append('    Process: frxst_pts was empty. Constant value raster created.')
        arcpy.AddMessage(loglines[-1])

        loglines = makeoutncfile(arcpy, constraster, 'gw_basns.nc', 'basn_msk', projdir, loglines)
        loglines.append('    Process: gw_basns was empty. Constant value raster created.')
        arcpy.AddMessage(loglines[-1])

    # Process: Alter Channelgrid for reservoirs if necessary
    if in_lakes != '':
        arcpy, channelgrid, outRaster, loglines = add_reservoirs(arcpy, channelgrid, in_lakes, flac, projdir, fill2, cellsize2, sr2, loglines)
        loglines = makeoutncfile(arcpy, outRaster, 'LAKEGRID.nc', 'LAKEGRID', projdir, loglines)
    else:
        loglines = makeoutncfile(arcpy, constraster, 'LAKEGRID.nc', 'LAKEGRID', projdir, loglines)

    # Process: Output Channelgrid
    loglines = makeoutncfile(arcpy, channelgrid, 'CHANNELGRID.nc', 'CHANNELGRID', projdir, loglines)

    del constraster, fdir, flac, strm, channelgrid, fill2
    # Process: Resample LU_INDEX grid to a higher resolution
    LU_INDEX2 = os.path.join(projdir, "LU_INDEX")
    arcpy.Resample_management(LU_INDEX, LU_INDEX2, cellsize2, "NEAREST")
    loglines = makeoutncfile(arcpy, LU_INDEX2, 'landuse.nc', 'landuse', projdir, loglines)
    
    # zip the folder
    nclist = ['longitude.nc',
                'latitude.nc',
                'topography.nc',
                'topography_org.nc',
                'flowdirection.nc',
                'flowacc.nc',
                'LAKEGRID.nc',
                'retdeprtfac.nc',
                'ovroughrtfac.nc',
                'str_order.nc',
                'frxst_pts.nc',
                'frxst_basns.nc',
                'gw_basns.nc',
                'CHANNELGRID.nc',
                'landuse.nc',
                'LINKID.nc',
                'gw_basns_geogrid.txt',
                'Route_Link.csv',
                'LAKEPARM.TBL',
                'streams.shp', 'streams.shx', 'streams.shp.xml', 'streams.sbx', 'streams.sbn', 'streams.prj', 'streams.dbf']
    zipper = zipUpFolder(arcpy, zipfile, projdir, out_zip, nclist)
    arcpy.Delete_management('mosaicprj')
    arcpy.Delete_management(mosprj)
    loglines.append('    Step 4 completed without error.')
    arcpy.AddMessage(loglines[-1])
    return loglines

def main():
    pass

if __name__ == '__main__':
    main()
