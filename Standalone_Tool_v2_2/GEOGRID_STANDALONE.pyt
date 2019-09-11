# *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
# ** Copyright UCAR (c) 1992 - 2012
# ** University Corporation for Atmospheric Research(UCAR)
# ** National Center for Atmospheric Research(NCAR)
# ** Research Applications Laboratory(RAL)
# ** P.O.Box 3000, Boulder, Colorado, 80307-3000, USA
# ** 2012/1/5 10:12:13
# *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

import arcpy
import os, sys, math
import shutil
sys.dont_write_bytecode = True
import wrf_hydro_functions
import time

# Find out if we are in 32-bit or 64-bit
if sys.maxsize > 2**32:
    bit64 = True
else:
    bit64 = False

# Turn functions script and TauDEM script into Project data so it gets copied up to the server directory
wrfhydro_script = r'E:\Projects\Gochis\TOOLS\GEOGRID_PREPROCESSOR_CURRENT\wrf_hydro_functions.py'

class Toolbox(object):
    def __init__(self):
        """Define the toolbox (the name of the toolbox is the name of the
        .pyt file)."""
        self.label = "GEOGRID_STANDALONE"
        self.alias = ""
        self.description = "This is a standalone processing tool for WRF-HYDRO."

        # List of tool classes associated with this toolbox
        self.tools = [ProcessGeogridFile, ExportGrid, ExamineOutputs, ExportPRJ, GenerateLatLon, DomainShapefile]

class ProcessGeogridFile(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "ProcessGeogridFile"
        self.description = "This tool takes an input WRF Geogrid file in NetCDF format" + \
                           " and uses the HGT_M grid and an input high-resolution elevation grid" + \
                           "to produce a high-resolution hydrologically processed output."
        #self.canRunInBackground = False
        self.canRunInBackground = True
        self.category = "Processing"

    def getParameterInfo(self):
        """Define parameter definitions"""

        in_nc = arcpy.Parameter(
            displayName="Input GEOGRID File",
            name="in_nc",
            datatype="File",
            parameterType="Required",
            direction="Input")

        in_csv = arcpy.Parameter(
            displayName="Forecast Points (CSV)",
            name="in_csv",
            datatype="File",
            parameterType="Optional",
            direction="Input")
        # To define a file filter that includes .csv and .txt extensions,
        #  set the filter list to a list of file extension names
        in_csv.filter.list = ['csv']

        basin_mask = arcpy.Parameter(
            displayName="Mask CHANNELGRID to basins?",
            name="basin_mask",
            datatype="Boolean",
            parameterType="Optional",
            direction="Input")
        basin_mask.value = False
        # datatype="GPBoolean" keyword can be used at 10.1 SP1

        RB_routing = arcpy.Parameter(
            displayName="Create reach-based routing files?",
            name="RB_routing",
            datatype="Boolean",
            parameterType="Optional",
            direction="Input")
        RB_routing.value = False
        # datatype="GPBoolean" keyword can be used at 10.1 SP1

        Lake_routing = arcpy.Parameter(
            displayName="Create lake parameter file?",
            name="Lake_routing",
            datatype="Boolean",
            parameterType="Optional",
            direction="Input")
        Lake_routing.value = False
        # datatype="GPBoolean" keyword can be used at 10.1 SP1

        in_reservoirs = arcpy.Parameter(
            displayName="Reservoirs Shapefile or Feature Class",
            name="in_reservoirs",
            datatype="DEFeatureClass",
            parameterType="Optional",
            direction="Input")

        in_raster = arcpy.Parameter(
            displayName="Input Raster",
            name="in_raster",
            datatype="Raster Dataset",
            parameterType="Required",
            direction="Input")
        # datatype="DERasterDataset" keyword can be used at 10.1 SP1

        cellsize = arcpy.Parameter(
            displayName="Regridding Factor",
            name="cellsize",
            datatype="Long",
            parameterType="Required",
            direction="Input")
        cellsize.value = 15
        # datatype="GPLong" keyword can be used at 10.1 SP1

        threshold = arcpy.Parameter(
            displayName="Number of pixels to define stream",
            name="threshold",
            datatype="Long",
            parameterType="Required",
            direction="Input")
        threshold.value = 50
        # datatype="GPLong" keyword can be used at 10.1 SP1

        ovroughrtfac_val = arcpy.Parameter(
            displayName="OVROUGHRTFAC Value",
            name="ovroughrtfac_val",
            datatype="Any Value",
            parameterType="Required",
            direction="Input")
        ovroughrtfac_val.value = '0.3'
        ovroughrtfac_val.category = "Parameter Values"

        retdeprtfac_val = arcpy.Parameter(
            displayName="RETDEPRTFAC Value",
            name="retdeprtfac_val",
            datatype="Any Value",
            parameterType="Required",
            direction="Input")
        retdeprtfac_val.value = '0.0'
        retdeprtfac_val.category = "Parameter Values"

        out_zip = arcpy.Parameter(
            displayName="Output ZIP File",
            name="out_zip",
            datatype="File",
            parameterType="Required",
            direction="Output")
        out_zip.value = 'WRF_Hydro_routing_grids.zip'

        parameters = [in_nc, in_csv, basin_mask, RB_routing, Lake_routing, in_reservoirs, in_raster, cellsize, threshold, ovroughrtfac_val, retdeprtfac_val, out_zip]
        return parameters

    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True

    def isLicensed(self):
        """Allow the tool to execute, only if the ArcGIS Spatial Analyst extension
        is available."""
        try:
            if arcpy.CheckExtension("Spatial") != "Available":
                raise Exception
        except Exception:
            return False                                                            # tool cannot be executed
        return True                                                                 # tool can be executed

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""

        # Only activate masking if a CSV file has been input
        if not parameters[1].altered:
          parameters[2].enabled = False
        else:
          parameters[2].enabled = True

        # Only activate Lake input parameter if requesting lake CSV file
        if parameters[4].value == True:
            parameters[5].enabled = True
        else:
            parameters[5].enabled = False

        # Get the directory name from input parameter 0 and use it for the output
        if parameters[0].altered:                                               #and not parameters[0].hasBeenValidated
            workspace = os.path.dirname(parameters[0].value.value)
            parameters[11].value = os.path.join(workspace, parameters[11].value.value)
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""

        # Gather all necessary parameters
        in_nc = parameters[0].valueAsText
        in_csv = parameters[1].valueAsText
        basin_mask = parameters[2].valueAsText
        routing = parameters[3].valueAsText
        Lake_routing = parameters[4].valueAsText
        in_reservoir = parameters[5].valueAsText
        in_raster = parameters[6].valueAsText
        cellsize = parameters[7].valueAsText
        threshold = parameters[8].valueAsText
        ovroughrtfac_val = parameters[9].valueAsText
        retdeprtfac_val = parameters[10].valueAsText
        out_zip = parameters[11].valueAsText

        # Prepare output log file
        outtable = open(os.path.join(os.path.dirname(out_zip), os.path.basename(out_zip) + '.log'), "w")
        loglines = ['Begining processing on %s' %time.ctime()]
        loglines.append('64-bit background geoprocessing: %s' %bit64)
        tic = time.time()
        loglines.append('Input parameters:')
        for param in parameters:
            loglines.append('    Parameter: %s: %s' %(param.displayName, param.valueAsText))
        outtable.writelines("\n".join(loglines) + "\n")

        # Create scratch directory for temporary outputs
        projdir = os.path.dirname(out_zip) + os.sep + 'scratchdir'              # This is the only instance where we need the 'os' module!
        if not os.path.exists(projdir):
            os.makedirs(projdir)
        arcpy.env.overwriteOutput = True
        arcpy.env.workspace = projdir
        arcpy.env.scratchWorkspace = projdir

        # Set the input units to meters
        inunits = 'm'

        # Interpret the input for reservoir routing
        if Lake_routing == 'false':
            in_lakes = ''
        else:
            in_lakes = in_reservoir

        # Step 1 - Georeference geogrid file

        LU_INDEX, sr2, Projection_String, loglines = wrf_hydro_functions.georeference_geogrid_file(arcpy, in_nc, 'LU_INDEX')            # Process: Generate LU Index grid
        outtable.writelines("\n".join(loglines) + "\n")
        
        hgt_m_raster, sr2, Projection_String, loglines = wrf_hydro_functions.georeference_geogrid_file(arcpy, in_nc, 'HGT_M')
        outtable.writelines("\n".join(loglines) + "\n")


        # Step 2 - Create high resolution topography layers
        mosprj, cellsize1, cellsize2, loglines = wrf_hydro_functions.create_high_res_topogaphy(arcpy, in_raster, hgt_m_raster, cellsize, sr2, projdir)
        outtable.writelines("\n".join(loglines) + "\n")

        # Step 3 - Create latitude and longitude rasters
        loglines = wrf_hydro_functions.create_lat_lon_rasters(arcpy, projdir, mosprj)
        outtable.writelines("\n".join(loglines) + "\n")

        # Step 4 - Hyrdo processing functions
        loglines = wrf_hydro_functions.sa_functions(arcpy, basin_mask, mosprj, ovroughrtfac_val, retdeprtfac_val, projdir, in_csv, out_zip, threshold, inunits, LU_INDEX, cellsize1, cellsize2, routing, in_lakes) # , mosprj2,
        outtable.writelines("\n".join(loglines) + "\n")

        # Clean up and give finishing message
        #del LU_INDEX, hgt_m_raster
        #shutil.rmtree(projdir)
        loglines = ['Completed without error in %s seconds.\n' %(time.time()-tic)]
        arcpy.AddMessage(loglines[-1])
        outtable.write(loglines[-1])
        outtable.close()
        return

class ExportGrid(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "ExportGrid"
        self.description = "This tool takes an input WRF Geogrid file in NetCDF format" + \
                           " and uses the specified variable's grid to produce a raster."
        self.canRunInBackground = True
        self.category = "Utilities"

    def getParameterInfo(self):
        """Define parameter definitions"""

        # Input parameter
        in_nc = arcpy.Parameter(
            displayName="Input GEOGRID File",
            name="in_nc",
            datatype="File",
            parameterType="Required",
            direction="Input")

        # Fourth parameter
        var_name = arcpy.Parameter(
            displayName="Variable Name",
            name="var_name",
            datatype="GPString",
            parameterType="Required",
            direction="Input")
        var_name.filter.type = "ValueList"

        # Third parameter
        out_raster = arcpy.Parameter(
            displayName="Output Raster",
            name="out_raster",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Output")

        parameters = [in_nc, var_name, out_raster]
        return parameters

    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        if parameters[0].altered:
            in_nc_file = parameters[0].valueAsText

            # Establish an object for reading the input NetCDF file
            ncFP = arcpy.NetCDFFileProperties(in_nc_file)

            # Loop through global variables in NetCDF file to gather projection information
            ncVarNames = ncFP.getVariablesByDimension('west_east')
            ncMassgridNames = []
            for x in ncVarNames:
                mgridvar = ncFP.getAttributeValue(x, 'stagger')                 # Only use variables on Massgrid for now ('M')
                if mgridvar == 'M':
                    ncMassgridNames.append(x)
            parameters[1].filter.list = ncMassgridNames
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""

        # Set environments
        arcpy.env.overwriteOutput = True

        # Gather all necessary parameters
        in_nc = parameters[0].valueAsText
        Variable = parameters[1].valueAsText
        out_raster = parameters[2].valueAsText

        # Use wrf_hydro_functions to perform process
        nc_raster, sr2, Projection_String, loglines = wrf_hydro_functions.georeference_geogrid_file(arcpy, in_nc, Variable)

        # Set environments and save
        arcpy.env.outputCoordinateSystem = sr2
        nc_raster.save(out_raster)
        arcpy.DefineProjection_management(out_raster, sr2)
        del nc_raster

        arcpy.AddMessage('    Process completed without error.')
        arcpy.AddMessage('    Output Raster: %s' %out_raster)
        return

class ExamineOutputs(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "ExamineOutputs"
        self.description = "This tool takes the output zip file from the ProcessGeogrid script" + \
                           "and creates a raster from each output NetCDF file." + \
                           "" + \
                           "The Input should be a .zip file that was created using the WRF Hydro pre-" + \
                           "processing tools.  The Output Folder parameter should be set to a non-existent " +\
                           "folder location.  The tool will create the folder which will contain the results."
        self.canRunInBackground = True
        self.category = "Utilities"

    def getParameterInfo(self):
        """Define parameter definitions"""

        # Input parameter
        in_zip = arcpy.Parameter(
            displayName="Input ZIP File",
            name="in_zip",
            datatype="File",
            parameterType="Required",
            direction="Input")

        # Output parameter
        out_folder = arcpy.Parameter(
            displayName="Output Folder",
            name="out_folder",
            datatype="DEFolder",
            parameterType="Required",
            direction="Output")
        out_folder.defaultEnvironmentName = "workspace"

        parameters = [in_zip, out_folder]
        return parameters

    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""

        # Initiate input and output parameters
        in_zip = parameters[0].valueAsText
        out_folder = parameters[1].valueAsText

        # Use wrf_hydro_functions to perform process
        wrf_hydro_functions.Examine_Outputs(arcpy, in_zip, out_folder)
        return

class ExportPRJ(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "ExportPRJ"
        self.description = "This tool takes an input WRF Geogrid file in NetCDF format" + \
                           " and uses the specified variable's projection parameters" + \
                           " to produce a projection file."
        self.canRunInBackground = True
        self.category = "Utilities"

    def getParameterInfo(self):
        """Define parameter definitions"""

        # Input parameter
        in_nc = arcpy.Parameter(
            displayName="Input GEOGRID File",
            name="in_nc",
            datatype="File",
            parameterType="Required",
            direction="Input")

        # Second parameter
        out_prj = arcpy.Parameter(
            displayName="Output Projection File (.prj)",
            name="out_prj",
            datatype="DEPrjFile",
            parameterType="Required",
            direction="Output")

        parameters = [in_nc, out_prj]
        return parameters

    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""

        ##if parameters[0].altered:
        ##    in_nc_file = parameters[0].valueAsText
        ##
        ##    # Set output workspace and default name
        ##    workspace = os.path.dirname(in_nc_file)
        ##    parameters[1].value = os.path.join(workspace, in_nc_file + '.prj')

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""

        # Set environments
        arcpy.env.overwriteOutput = True

        # Gather all necessary parameters
        in_nc = parameters[0].valueAsText
        out_prj = parameters[1].valueAsText
        Variable = 'HGT_M'

        # Use wrf_hydro_functions to perform process
        nc_raster, sr2, Projection_String, loglines = wrf_hydro_functions.georeference_geogrid_file(arcpy, in_nc, Variable)
        del nc_raster, sr2

        # Optional save .prj file
        with open(out_prj, 'w') as prj_file:
             prj_file.write(Projection_String)

        arcpy.AddMessage('    Process completed without error.')
        arcpy.AddMessage('    Output Projection File: %s' %out_prj)
        return

class GenerateLatLon(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "GenerateLatLon"
        self.description = "This tool takes an input raster (most likely produced" + \
                           " using the ExportGrid tool) and uses that grid to produce" + \
                           " latitude and longitude netcdf file and ESRI GRID outputs."
        self.canRunInBackground = True
        self.category = "Utilities"

    def getParameterInfo(self):
        """Define parameter definitions"""

        # Input parameter
        in_raster = arcpy.Parameter(
            displayName="Input Raster",
            name="in_raster",
            datatype="Raster Dataset",
            parameterType="Required",
            direction="Input")
        # datatype="DERasterDataset" keyword can be used at 10.1 SP1

        # Output parameter
        out_folder = arcpy.Parameter(
            displayName="Output Folder",
            name="out_folder",
            datatype="DEFolder",
            parameterType="Required",
            direction="Input")
        out_folder.defaultEnvironmentName = "workspace"

        parameters = [in_raster, out_folder]
        return parameters

    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""

        # Set environments
        arcpy.env.overwriteOutput = True

        # Gather all necessary parameters
        inraster = parameters[0].valueAsText
        projdir = parameters[1].valueAsText
        arcpy.AddMessage('Input Raster Dataset: %s' %inraster)
        arcpy.AddMessage('Directory to be used for outputs: %s' %projdir)

        # Step 3 - Create latitude and longitude rasters
        in_raster = arcpy.Raster(inraster)
        loglines = wrf_hydro_functions.create_lat_lon_rasters(arcpy, projdir, in_raster)

        arcpy.AddMessage('    Process completed without error.')
        return

class DomainShapefile(object):

    """This function taks the WRF GEOGRID file and outputs a polygon shapefile using the bounding coordinate information."""

    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "DomainShapefile"
        self.description = "This tool takes an WRF Geogrid file and creates a single" + \
                           " polygon shapefile that makes up the boundary of the domain" + \
                           " of the M-grid (HGT_M, for example)."
        self.canRunInBackground = True
        self.category = "Utilities"

    def getParameterInfo(self):
        """Define parameter definitions"""

        # Input parameter
        in_nc = arcpy.Parameter(
            displayName="Input GEOGRID File",
            name="in_nc",
            datatype="File",
            parameterType="Required",
            direction="Input")

        # Output parameter
        out_shp = arcpy.Parameter(
            displayName="Output Shapefile",
            name="out_shp",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Output")

        parameters = [in_nc, out_shp]
        return parameters

    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        try:
            if arcpy.CheckExtension("Spatial") != "Available":
                raise Exception
        except Exception:
            return False                                                            # tool cannot be executed
        return True                                                                 # tool can be executed

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""

        # Gather all necessary parameters
        in_nc = parameters[0].valueAsText
        out_shp = parameters[1].valueAsText

        # Create scratch directory for temporary outputs
        projdir = os.path.dirname(in_nc)
        arcpy.env.overwriteOutput = True
        arcpy.env.workspace = projdir
        arcpy.env.scratchWorkspace = projdir

        # Get coarse grid raster
        hgt_m_raster, sr2, Projection_String, loglines = wrf_hydro_functions.georeference_geogrid_file(arcpy, in_nc, 'HGT_M')
        del Projection_String, loglines

        # Use coarse grid raster to create shapefile
        wrf_hydro_functions.domain_shapefile(arcpy, hgt_m_raster, out_shp, sr2)

        # Delete intermediate files in scratch dir
        del hgt_m_raster, sr2

        return parameters
