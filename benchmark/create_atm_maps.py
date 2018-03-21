#!/usr/bin/python
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np

output_dir = 'atm_maps'
#sub_dirs = ['temperature','v_velocity','w_velocity','velocity_mag', 'water_vapour', 'precipitation', 'precipitable_water']

def create_maps(directory):
    index = 6
    if directory == 'temperature':
        index = 6
    elif directory == 'v_velocity':
        index = 3
    elif directory == 'w_velocity':
        index = 4
    elif directory == 'velocity_mag':
        index = 5
    elif directory == 'water_vapour':
        index = 7
    elif directory == 'precipitation':
        index = 8
    elif directory == 'precipitable_water':
        index = 9
    elif  directory == 'topography':
        index = 10
    elif  directory == 'topo_simon':
        index = 11

    for time in range(0,150,10):
        if index < 10:
            data = np.genfromtxt('./output/[{}Ma_Golonka.xyz]_PlotData_Atm.xyz'.format(time),skip_header=1)
            for d in data:
                d[1]=90-d[1]
            x = data[:,0]
            y = data[:,1]
            z = data[:,index]
        elif index == 10:
            data = np.genfromtxt('../data/Golonka_Smoothed/{}Ma_Golonka.xyz'.format(time))
            x = data[:,0]
            y = data[:,1]
            z = data[:,2]
        
        elif index == 11:
            data = np.genfromtxt('./paleotopo_grids/{}Ma.xyz'.format(time))
            x = data[:,0]
            y = data[:,1]
            z = data[:,2]

        topo = data[:,2]

        plt.figure(figsize=(15, 8))

        m = Basemap(llcrnrlon=-180,llcrnrlat=-90,urcrnrlon=180,urcrnrlat=90,projection='kav7', lon_0=0)  

        if index != 10 and index != 11:
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

        if index != 10 and index != 11:
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
    #output_dir = 'atm_maps'
    #v-velocity(m/s), w-velocity(m/s), velocity-mag(m/s), temperature(Celsius), water_vapour(g/kg), precipitation(mm), precipitable water(mm)
    #sub_dirs = ['temperature','v_velocity','w_velocity','velocity_mag', 'water_vapour', 'precipitation', 'precipitable_water', 'topography', 'topo_simon']
    sub_dirs = ['topo_simon']
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for d in sub_dirs:
        if not os.path.exists(output_dir+'/'+d):
            os.makedirs(output_dir+'/'+d)

        create_maps(d)




