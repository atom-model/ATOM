#!/usr/bin/python
# -*- coding: utf-8 -*-
import os, sys, shutil
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np

from draw_atm_velocities  import draw_velocity
from draw_topo import draw_topography

map_cfg = {
    'temperature': (6, -20, 30, u'Temperature (Celsius °C)'),
    'v_velocity': (3, -7, 7, 'v_velocity (m/s)'),
    'w_velocity': (4, -8, 9, 'w_velocity (m/s)'),
    'water_vapour': (7, 0, 20, 'Water Vapour (g/kg)'),
    'precipitation': (8, 0, 800, 'Precipitation (mm/yr)'),
    'precipitable_water': (9, 0, 30, 'Precipitable Water (mm)'),
    'evaporation' : (10, 0, 7, 'Evaporation (mm/d)')
}

def create_maps(directory, start_time, end_time, time_step, output_dir, data_dir, topo_dir, topo_suffix):
    for time in range(start_time, end_time + 1, time_step):
        if 'topography' == directory:
            draw_topography(time, output_dir + '/topography/', topo_dir, topo_suffix)
            continue

        if 'velocity' == directory:
            draw_velocity(time, output_dir + "/velocity/", data_dir)
            continue

        if directory not in map_cfg:
            print('unrecognized component: ' + directory)
            return

        data = np.genfromtxt(data_dir + '/[{0}Ma_{1}.xyz]_PlotData_Atm.xyz'.format(time, topo_suffix), skip_header=1)
        index = map_cfg[directory][0]
        x = data[:,0]
        y = data[:,1]
        z = data[:,index]
        if directory == 'precipitation':
            z=z*365

        topo = data[:,2]

        plt.figure(figsize=(15, 8))

        m = Basemap(llcrnrlon=-180,llcrnrlat=-90,urcrnrlon=180,urcrnrlat=90,projection='kav7', lon_0=0)  

        xx = x
        x, topo = m.shiftdata(xx, datain = topo, lon_0=0)
        x, z = m.shiftdata(xx, datain = z, lon_0=0)

        xi, yi = m(x, y)
        v_min = map_cfg[directory][1]
        v_max = map_cfg[directory][2] 

        #cs = m.scatter(xi, yi, marker='.', c=z, alpha=0.5, lw=0, vmin=v_min, vmax=v_max)

        zz = np.zeros((181, 361))
        for i in range(181):
            for j in range(361):
                zz[180-i][j] = z[j*181+i]
        img_data = m.transform_scalar(zz, np.arange(-180,180),np.arange(-90,90),361,181)
        cs = m.imshow(img_data,alpha=0.5, vmin=v_min, vmax=v_max, cmap='jet')

        m.contour( xi.reshape((361,181)), yi.reshape((361,181)), topo.reshape((361,181)),
                            colors ='k', linewidths= 0.3 )

        m.drawparallels(np.arange(-90., 90., 10.), labels=[1,0,0,0], fontsize=10)
        m.drawmeridians(np.arange(-180., 180., 45.), labels=[0,0,0,1], fontsize=10)
        #m.drawcoastlines()   

        cbar = m.colorbar(cs, location='bottom', pad="10%", label=map_cfg[directory][3])
        plt.title("{1} at {0}Ma".format(time, directory.capitalize()))
        plt.savefig(output_dir+'/'+directory+'/{0}_Ma_{1}.png'.format(time, directory), bbox_inches='tight')
        print(output_dir+'/'+directory+'/{0}_Ma_{1}.png has been saved!'.format(time, directory))
        #plt.show()
        plt.close()


def create_all_maps(sub_dirs, start_time, end_time, time_step, output_dir, data_dir, topo_dir, topo_suffix):
    answer = ""
    if os.path.exists(output_dir):
        while False:
            answer = raw_input("Do you want to delete the existing output folder {} [Y/N]? ".format(
                os.path.abspath(output_dir))).lower()
            if answer in ["y", "n"]:
                break
        answer = "y"
        if answer == "y":
            shutil.rmtree(output_dir)
            print("Delete the existing " + os.path.abspath(output_dir))
            os.makedirs(output_dir)
            print("Create a new " + os.path.abspath(output_dir))
    else:
        os.makedirs(output_dir)
        print("Create " + os.path.abspath(output_dir))

    print('')
    
    try:
        for d in sub_dirs:
            if not os.path.exists(output_dir+'/'+d):
                os.makedirs(output_dir+'/'+d)
            create_maps(d, start_time, end_time, time_step, output_dir, data_dir, topo_dir, topo_suffix)
    except:
        import traceback
        traceback.print_exc()

    print("The data being used to create the maps are in " + os.path.abspath(data_dir))
    print("The paleo-topography dir is: " + os.path.abspath(topo_dir))
    print("The paleo-topography file suffix is: " + topo_suffix)
    print("The start time is: " + str(start_time) )
    print("The end time is: " + str(end_time) )
    print("The time step is: " + str(time_step) )
    print("The output dir is: " + os.path.abspath(output_dir) )
    print("Creating atm maps complete!")

if  __name__ == "__main__":

    data_dir = '../benchmark/output'
    topo_dir = '../data/Paleotopography_bathymetry/Golonka_rev210/'
    topo_suffix = 'Golonka'
    start_time = 0
    end_time = 10
    time_step = 5
    output_dir = './atm_maps'

    # v-velocity(m/s), w-velocity(m/s), velocity-mag(m/s), temperature(Celsius), water_vapour(g/kg), 
    # precipitation(mm), precipitable water(mm)
    sub_dirs = ['temperature','v_velocity','w_velocity', 'water_vapour', 
        'precipitation', 'precipitable_water', 'topography', 'velocity', 'evaporation']

    try:
        start_time = int(sys.argv[1])
        end_time  = int(sys.argv[2])
        time_step = int(sys.argv[3])
        data_dir = sys.argv[4]
        topo_suffix = sys.argv[5]
        #topo_dir = sys.argv[6]
    except:
        pass

    create_all_maps(sub_dirs, start_time, end_time, time_step, output_dir, data_dir, topo_dir, topo_suffix)

    print("Use './create_atm_maps.py 0 100 5 ../benchmark/output Ma_smooth' to override the default settings")
