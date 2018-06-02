#!/usr/bin/python
# -*- coding: utf-8 -*-
import os, sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np

data_dir = './output'
start_time = 0
end_time = 100
time_step = 5

output_dir = 'hyd_maps'
#'salinity', 'bottom_water', 'upwelling', 'downwelling
def create_maps(directory):
    index = 6
    if directory == 'temperature':
        index = 6
    elif directory == 'v_velocity':
        index = 3
    elif directory == 'w_velocity':
        index = 4
    elif directory == 'salinity':
        index = 7
    elif directory == 'bottom_water':
        index = 8
    elif directory == 'upwelling':
        index = 9
    elif directory == 'downwelling':
        index = 10


    for time in range(start_time, end_time, time_step):
        if index < 11:
            data = np.genfromtxt(data_dir + '/[{}Ma_Golonka.xyz]_PlotData_Hyd.xyz'.format(time),skip_header=1)
            for d in data:
                d[1]=90-d[1]
            x = data[:,0]
            y = data[:,1]
            z = data[:,index]

        topo = data[:,2]

        plt.figure(figsize=(15, 8))

        m = Basemap(llcrnrlon=-180,llcrnrlat=-90,urcrnrlon=180,urcrnrlat=90,projection='kav7', lon_0=0)  

        xx = x
        x, topo = m.shiftdata(xx, datain = topo, lon_0=0)
        x, z = m.shiftdata(xx, datain = z, lon_0=0)

        xi, yi = m(x, y)
        v_min=None
        v_max=None
        if index == 6: #tempetature
            v_min = -60
            v_max = 35
        elif index == 8:
            v_min = 0
            v_max = 7

        cs = m.scatter(xi, yi, marker='.', c=z, alpha=0.5, lw=0, vmin=v_min, vmax=v_max)

        m.contour( xi.reshape((361,181)), yi.reshape((361,181)), topo.reshape((361,181)),
                            colors ='k', linewidths= 0.3 )

        m.drawparallels(np.arange(-90., 90., 10.), labels=[1,0,0,0], fontsize=10)
        m.drawmeridians(np.arange(-180., 180., 45.), labels=[0,0,0,1], fontsize=10)
        #m.drawcoastlines()   

        cbar = m.colorbar(cs, location='bottom', pad="10%")
        plt.title("{0}Ma {1}".format(time,directory))
        plt.savefig(output_dir+'/'+directory+'/{0}Ma_{1}.png'.format(time, directory), bbox_inches='tight')
        print '{0}Ma_{1}.png has been saved!'.format(time, directory)
        #plt.show()
        plt.close()

if  __name__ == "__main__":
    sub_dirs = ['temperature','v_velocity','w_velocity', 'salinity', 'bottom_water', 'upwelling', 'downwelling']

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    try:
        start_time = int(sys.argv[1])
        end_time  = int(sys.argv[2])
        time_step = int(sys.argv[3])
        data_dir = sys.argv[4]
        for d in sub_dirs:
            if not os.path.exists(output_dir+'/'+d):
                os.makedirs(output_dir+'/'+d)

            create_maps(d)
    except:
        print("Usage: ./create_hyd_maps.py 0 100 5")
        import traceback
        traceback.print_exc()






