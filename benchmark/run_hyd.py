#!/usr/bin/env python

import sys, random, argparse
sys.path.append('../reconstruction')
sys.path.append('../utils')

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np

from reconstruct_atom_data import *
from pyatom import Atmosphere, Hydrosphere
import create_atm_maps, create_hyd_maps

hyd_model = Hydrosphere()

hyd_model.load_config( './config_hyd.xml' )

start_time = hyd_model.time_start
end_time = hyd_model.time_end
time_step = hyd_model.time_step
print(('{0}  {1}  {2}'.format(start_time,end_time, time_step)))
times = list(range(start_time, end_time+1, time_step))

atom_output_dir = hyd_model.output_path

BATHYMETRY_SUFFIX = 'Ma_smooth.xyz'

for t in range(len(times)):
    time = times[t]
#    atm_model.run_time_slice(time)
    hyd_model.run_time_slice(time)
    #if False:
#    if t<len(times)-1:
#        reconstruct_temperature(time,times[t+1], BATHYMETRY_SUFFIX) 
#        reconstruct_precipitation(time,times[t+1], BATHYMETRY_SUFFIX)
#        reconstruct_salinity(time,times[t+1], BATHYMETRY_SUFFIX)
#        reconstruct_wind_v(time,times[t+1], BATHYMETRY_SUFFIX)
#        reconstruct_wind_w(time,times[t+1], BATHYMETRY_SUFFIX)
try:
    topo_dir = '../data/topo_grids/'
    topo_suffix = 'smooth'
    
#    atm_map_output_dir = './atm_maps'
    hyd_map_output_dir = './hyd_maps'

    # v-velocity(m/s), w-velocity(m/s), velocity-mag(m/s), temperature(Celsius), water_vapour(g/kg), 
    # precipitation(mm), precipitable water(mm)
#    atm_sub_dirs = ['temperature','v_velocity','w_velocity', 'water_vapour', 'precipitation', 
#                'precipitable_water', 'topography', 'velocity', 'evaporation']

#    create_atm_maps.create_all_maps(atm_sub_dirs, start_time, end_time, time_step, atm_map_output_dir, 
#            atom_output_dir, topo_dir, topo_suffix)

    hyd_sub_dirs = ['temperature','v_velocity','w_velocity', 'salinity', 'Ekman_pumping', 
            'upwelling', 'downwelling', 'velocity']
    
    create_hyd_maps.create_all_maps(hyd_sub_dirs, start_time, end_time, time_step, hyd_map_output_dir,
            atom_output_dir, topo_dir, topo_suffix)
except:
    import traceback
    traceback.print_exc() 
