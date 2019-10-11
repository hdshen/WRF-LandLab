# Written by Brigid Lynch (Indiana University, Dept. of Earth and
# Atmospheric Sciences. This is an example script showing how to build a simple rectangular
# mountain with 2 flanks and channel networks.

from landlab.components.flow_routing.route_flow_dn import FlowRouter
from landlab.components.stream_power.stream_power import StreamPowerEroder
from landlab.components.diffusion.diffusion import LinearDiffuser
from landlab import RasterModelGrid, CLOSED_BOUNDARY, FIXED_VALUE_BOUNDARY
from landlab import ModelParameterDictionary
from landlab.io.netcdf import write_netcdf
import numpy as np

#get needed properties to build grid:
input_file = 'create_topo.txt'

#initialize an object that will supply the parameters:
inputs = ModelParameterDictionary(input_file)

#load grid properties
nrows = inputs.read_int('nrows') 
ncols = inputs.read_int('ncols') 
dx = inputs.read_float('dx') 
dt = inputs.read_float('dt') 
total_t = inputs.read_float('total_time') 
uplift_rate = inputs.read_float('uplift_rate') 

nt = int(total_t // dt) #this is how many loops we'll need
uplift_per_step = uplift_rate * dt

#build the grid
mg = RasterModelGrid((nrows, ncols), dx)
z = mg.add_zeros('node', 'topographic__elevation')
# add some roughness as the initial topography, this lets "natural" channel planforms arise
initial_roughness = np.random.rand(z.size)/100000.
z += initial_roughness

#set the boundary conditions - here we have closed boundaries at the top/bottom and fixed boundaries at the left/right
for edge in (mg.nodes_at_left_edge, mg.nodes_at_right_edge):
    mg.status_at_node[edge] = FIXED_VALUE_BOUNDARY
for edge in (mg.nodes_at_top_edge, mg.nodes_at_bottom_edge):
    mg.status_at_node[edge] = CLOSED_BOUNDARY

#instantiate the components
fr = FlowRouter(mg, input_file)
sp = StreamPowerEroder(mg, input_file)
lin_diffuse = LinearDiffuser(mg, input_file)

#run the model
elapsed_time = 0. #total time in simulation
for i in range(nt):
    lin_diffuse.run_one_step(dt) 
    fr.run_one_step() # route_flow isn't time sensitive, so it doesn't take dt as input
    sp.run_one_step(dt)
    mg.at_node['topographic__elevation'][mg.core_nodes] += uplift_per_step # add the uplift

    # if you want to see the evolution of the topogrpahy through time (say to check for steady state) you can output a topography at each time step
    #write_netcdf(('topography_output'+ str(elapsed_time) +'.nc'), mg, names='topographic__elevation')

    elapsed_time += dt

    if i % 100 == 0:
        print('completed loop %d' %i)

#save the resultant topography as a netcdf file
write_netcdf('landlab_topo_n.nc', mg, names='topographic__elevation')
