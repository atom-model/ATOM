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
from draw_temperature_plot import draw_temperature_plot
from draw_precipitation_map import draw_precipitation_map
from draw_precipitation_evaporation_diff_map import draw_precipitation_evaporation_diff_map

def main(maps_only=False):
    atm_model = Atmosphere()
    hyd_model = Hydrosphere()

    atm_model.load_config('./config_atm.xml')
    hyd_model.load_config('./config_hyd.xml')

    start_time = atm_model.time_start
    end_time = atm_model.time_end
    time_step = atm_model.time_step
    print('{0}  {1}  {2}'.format(start_time,end_time, time_step))
    times = range(start_time, end_time+1, time_step)

    atom_output_dir = atm_model.output_path

    BATHYMETRY_SUFFIX = 'Ma_smooth.xyz'

    if not maps_only:
        for t in range(len(times)):
            time = times[t]
            atm_model.run_time_slice(time)
            hyd_model.run_time_slice(time)
            if t<len(times)-1:
                reconstruct_temperature(time,times[t+1], BATHYMETRY_SUFFIX) 
                #reconstruct_precipitation(time,times[t+1], BATHYMETRY_SUFFIX)
                reconstruct_salinity(time,times[t+1], BATHYMETRY_SUFFIX)
                reconstruct_wind_v(time,times[t+1], BATHYMETRY_SUFFIX)
                reconstruct_wind_w(time,times[t+1], BATHYMETRY_SUFFIX)
    try:
        topo_dir = '../data/topo_grids/'
        topo_suffix = 'smooth'
    
        atm_map_output_dir = './output/atm_maps'
        hyd_map_output_dir = './output/hyd_maps'

        # v-velocity(m/s), w-velocity(m/s), velocity-mag(m/s), temperature(Celsius), water_vapour(g/kg), 
        # precipitation(mm), precipitable water(mm)
        atm_sub_dirs = ['temperature','v_velocity','w_velocity', 'water_vapour', 
                'precipitable_water', 'topography', 'velocity', 'evaporation']
#        atm_sub_dirs = ['temperature','v_velocity','w_velocity', 'water_vapour', 
#                'precipitation', 'precipitable_water', 'topography', 'velocity', 'evaporation']

        create_atm_maps.create_all_maps(atm_sub_dirs, start_time, end_time, time_step, atm_map_output_dir, 
            atom_output_dir, topo_dir, topo_suffix)

        hyd_sub_dirs = ['temperature','v_velocity','w_velocity', 'salinity', 'ekman_pumping', 
            'upwelling', 'downwelling', 'velocity']
    
        create_hyd_maps.create_all_maps(hyd_sub_dirs, start_time, end_time, time_step, hyd_map_output_dir,
            atom_output_dir, topo_dir, topo_suffix)

        if not os.path.isdir(atm_map_output_dir+'/precipitation/'):
            os.mkdir(atm_map_output_dir+'/precipitation/')

        if not os.path.isdir(atm_map_output_dir+'/precipitation_evaporation_diff/'):
            os.mkdir(atm_map_output_dir+'/precipitation_evaporation_diff/')

        for time in times:
            draw_precipitation_map(time, './output/', output_dir=atm_map_output_dir +'/precipitation/', topo_suffix=topo_suffix)
            draw_precipitation_evaporation_diff_map(time, './output/', 
                output_dir=atm_map_output_dir+'/precipitation_evaporation_diff/', topo_suffix=topo_suffix)

        if end_time >= 100:
            draw_temperature_plot(lon=180,data_dir='./output/', output_dir='./output/')

    except:
        import traceback
        traceback.print_exc() 

if __name__ == "__main__":
    __description__ = \
    """Some description words here
    """

    parser = argparse.ArgumentParser(
        description = __description__,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('--maps-only', action='store_true', help='do not run models. create maps only.')
    args = parser.parse_args()
    main(maps_only = args.maps_only)
    #print(args.maps_only)

