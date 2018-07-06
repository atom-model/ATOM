import sys, random
sys.path.append('../reconstruction')
sys.path.append('../utils')

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np

from reconstruct_atom_data import *
from pyatom import Model
import create_atm_maps, create_hyd_maps

model = Model()

start_time = 0
end_time = 30
time_step = 5
times = range(start_time, end_time, time_step)

atom_output_dir = './output'

for t in range(len(times)):
    time = times[t]
    model.run_atm( time, atom_output_dir, './config_atm.xml' )
    model.run_hyd( time, atom_output_dir, './config_hyd.xml' )
    #if False:
    if t<len(times)-1:
        reconstruct_temperature(time,times[t+1]) 
        reconstruct_precipitation(time,times[t+1])
        reconstruct_salinity(time,times[t+1])

try:
    topo_dir = '../data/Paleotopography_bathymetry/Golonka_rev210/'
    topo_suffix = 'Ma_Golonka'
    map_output_dir = './atm_maps'

    # v-velocity(m/s), w-velocity(m/s), velocity-mag(m/s), temperature(Celsius), water_vapour(g/kg), 
    # precipitation(mm), precipitable water(mm)
    sub_dirs = ['temperature','v_velocity','w_velocity', 'water_vapour', 'precipitation', 'precipitable_water', 'topography']

    create_atm_maps.create_all_maps(sub_dirs, start_time, end_time, time_step, map_output_dir, 
            atom_output_dir, topo_dir, topo_suffix)
except:
    import traceback
    traceback.print_exc() 
