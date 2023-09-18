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

import time
timestart_time = time.time()
from datetime import datetime

start_time_computer = datetime.now()
print('\n\n\n')
print('     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   Ocean General Circulation Model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \n') 
print('\n')
print('\n' + ' ... run_hyd.py: date and time at run time begin: \n' 
    + ' ... {0}'.format(start_time_computer))

hyd_model = Hydrosphere()

hyd_model.load_config( './config_hyd.xml' )

start_time = hyd_model.time_start
end_time = hyd_model.time_end
time_step = hyd_model.time_step

print('\n')
print('\n' + ' ... run_hyd.py: computed time range: \n' 
    + ' ... start time = {0} \n ... end time = {1} \n ... time step = {2}'
    .format(start_time, end_time, time_step))

times = range(start_time, end_time+1, time_step)

atom_output_dir = hyd_model.output_path

BATHYMETRY_SUFFIX = 'Ma_smooth.xyz'

for t in range(len(times)):
    time = times[t]
    hyd_model.run_time_slice(time)

    if t<len(times)-1:
        reconstruct_temperature(time,times[t+1], BATHYMETRY_SUFFIX) 
#        reconstruct_precipitation(time,times[t+1], BATHYMETRY_SUFFIX)
#        reconstruct_salinity(time,times[t+1], BATHYMETRY_SUFFIX)
#        reconstruct_wind_v(time,times[t+1], BATHYMETRY_SUFFIX)
#        reconstruct_wind_w(time,times[t+1], BATHYMETRY_SUFFIX)

end_time_computer = datetime.now()
print('\n')
print('\n' + ' ... run_hyd.py: date and time at run time end: \n' 
    + ' ... {0}'.format(end_time_computer))
print('\n' + ' ... run_hyd.py: duration of run time in hours, minutes and seconds:     {}'
    .format(end_time_computer - start_time_computer))
print('\n\n\n')

try:
    topo_dir = '../data/topo_grids/'
    topo_suffix = 'smooth'
    
    hyd_map_output_dir = './hyd_maps'
    hyd_sub_dirs = ['temperature','v_velocity','w_velocity', 'salinity', 'Ekman_pumping', 
                    'upwelling', 'downwelling', 'velocity']
    
    create_hyd_maps.create_all_maps(hyd_sub_dirs, start_time, end_time, time_step, hyd_map_output_dir,
            atom_output_dir, topo_dir, topo_suffix)

except:
    import traceback
    traceback.print_exc() 
