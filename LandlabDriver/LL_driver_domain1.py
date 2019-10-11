# Written by Brigid Lynch (Indiana University, Dept. of Earth and Atmospheric Sciences. Example
# Landlab code that uses discharge output from WRF-Hydro to drive erosion in landscape. This
# example uses discharge and topography data generated for Domain 1 in our study.

from landlab.components.flow_routing.route_flow_dn import FlowRouter
from landlab.components.stream_power.stream_power import StreamPowerEroder
from landlab.components.diffusion.diffusion import LinearDiffuser
from landlab import ModelParameterDictionary
from landlab import RasterModelGrid, CLOSED_BOUNDARY, FIXED_VALUE_BOUNDARY
from landlab.io.netcdf import write_netcdf
import numpy as np
from netCDF4 import Dataset
import time


#get needed properties to build grid:
input_file = 'LL_driver_params.txt'


#initialize an object that will supply the parameters:
inputs = ModelParameterDictionary(input_file)

#load grid properties
x = inputs.read_int('x') 
y = inputs.read_int('y') 
dx = inputs.read_float('dx') 
dt = inputs.read_float('dt') 
run_time = inputs.read_float('run_time') 
uplift = inputs.read_float('uplift_rate') 
initial_slope = inputs.read_float('initial_slope') 
nt = int(run_time//dt) #this is how many loops we'll need


#Build the grid
mg = RasterModelGrid(x, y, dx) #sp

# Load topography output rotated from WRF-Hydro and add to topography field
# Since the discharge files were created on the WRF-Hydro topography, we use that topography as our input. The original Landlab topography and the WRF-Hydro topography are almost identical, with only slight differences created by the rotation.
topo=Dataset('topography_wrf_hydro.nc')
topography=topo.variables['TOPOGRAPHY'][:]
topographic__elevation=np.asarray(topography, dtype=float).ravel()
mg.add_field('node', 'topographic__elevation', topographic__elevation) #create the field


#set boundary conditions -- these should match those used in the topography creation driver
for edge in (mg.nodes_at_left_edge, mg.nodes_at_right_edge):
    mg.status_at_node[edge] = CLOSED_BOUNDARY
for edge in (mg.nodes_at_top_edge, mg.nodes_at_bottom_edge):
    mg.status_at_node[edge] = FIXED_VALUE_BOUNDARY

#instantiate the components
fr = FlowRouter(mg)
sp = StreamPowerEroder(mg, input_file)
lin_diffuse = LinearDiffuser(mg, input_file)

print( 'Running ...' )
time_on = time.time()

#text for discharge file names
streamflow = "streamflow"
nc = ".nc"

elapsed_time = 0. #total time in simulation

#create list of all discharge files to be used, repeat as many times as necessary for run time.
#Example: here we have a total of 20 discharge files -- we need one file for every time step and we have 100 timesteps. 
# We will repeat the discharge file namelist 5 times Run Time = 100 yrs dt = 1yr 
rep_num = int(5) #depends on how many discharge files/time steps you have
discharge_filename=[streamflow + repr(j) + nc for k in range(rep_num) for j in range(1,21)]

#Note: we have included the code below that we use in this study: We use a total of 7301 discharge files.
# we need one file for every time step and we have 10,000 timesteps. Run Time = 1,000,000 yrs dt = 100 yrs
#Repeat the discharge file namelist 2 times (14,602).
#rep_num = int(2) #depends on how many discharge files/ time steps you have
#discharge_filename=[streamflow + repr(j) + nc for k in range(rep_num) for j in range(1,7302)]


#run the model: load a new discharge file at each time step, use this discharge file in the stream power erosion component
while elapsed_time < run_time:
    num = int(elapsed_time/dt)
    
    #load discharge to a new variable my_Q
    print (discharge_filename[num])
    f=Dataset(discharge_filename[num])
    discharge=f.variables['streamflow'][:]
    my_Q=np.asarray(discharge, dtype=float).ravel()
    
    #run diffusion component
    lin_diffuse.run_one_step(dt,deposit='false') #turn deposition off
    #run flow routing component
    fr.run_one_step()
    
    #add my_Q variable to landlab discharge
    _ = mg.add_field('node','surface_water__discharge', my_Q, noclobber=False)
    #run stream power component with my_Q
    sp.erode(mg, dt,Q_if_used=my_Q)
    
   
    #add uplift
    mg.at_node['topographic__elevation'][mg.core_nodes] +=uplift*dt
   
    #record topography through time (optional)
    #write_netcdf(('topography_output'+ str(elapsed_time) +'.nc'), mg, names='topographic__elevation')

    elapsed_time += dt
time_off = time.time()

#save the resultant topography as a netcdf file
write_netcdf('topography_final.nc', mg, names='topographic__elevation')


