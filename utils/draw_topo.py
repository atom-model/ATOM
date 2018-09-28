#!/usr/bin/python
# -*- coding: utf-8 -*-
import os, sys, shutil
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np

def draw_topography(time, output_dir, topo_dir, topo_suffix):
    data = np.genfromtxt(topo_dir + '/{0}Ma_{1}.xyz'.format(time, topo_suffix))
    x = data[:,0]
    y = data[:,1]
    z = data[:,2]
        
    plt.figure(figsize=(15, 8))

    m = Basemap(llcrnrlon=-180,llcrnrlat=-90,urcrnrlon=180,urcrnrlat=90,projection='kav7', lon_0=0)  

    xi, yi = m(x, y)
    v_min = -5000
    v_max = 5000

    #cs = m.scatter(xi, yi, marker='.', c=z, alpha=0.5, lw=0, vmin=v_min, vmax=v_max)

    zz = z.reshape((181, 361))
    zz = np.flipud(zz)
    #for i in range(181):
    #    for j in range(361):
    #        zz[180-i][j] = z[j*181+i]
    img_data = m.transform_scalar(zz, np.arange(-180,180),np.arange(-90,90),361,181)
    cs = m.imshow(img_data,alpha=0.5, vmin=v_min, vmax=v_max)


    m.drawparallels(np.arange(-90., 90., 10.), labels=[1,0,0,0], fontsize=10)
    m.drawmeridians(np.arange(-180., 180., 45.), labels=[0,0,0,1], fontsize=10)

    cbar = m.colorbar(cs, location='bottom', pad="10%", label='Topography (m)')
    
    plt.title("{1} at {0}Ma".format(time, 'Topography'))
    plt.savefig(output_dir+'/{0}Ma_{1}.png'.format(time, 'topography'), bbox_inches='tight')
    print(output_dir + '/{0}Ma_{1}.png has been saved!'.format(time, 'topography'))
    #plt.show()
    plt.close()

if  __name__ == "__main__":

    topo_dir = '../data/Paleotopography_bathymetry/Golonka_rev210/'
    topo_suffix = 'Golonka'
    start_time = 0
    end_time = 10
    time_step = 5
    output_dir = './topography'

    try:
        start_time = int(sys.argv[1])
        end_time  = int(sys.argv[2])
        time_step = int(sys.argv[3])
        topo_suffix = sys.argv[4]
        topo_dir = sys.argv[5]
    except:
        print("Usage: python " + sys.argv[0] +  ' 0 20 5 Golonka ../data/Paleotopography_bathymetry/Golonka_rev210/')

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    for time in range(start_time, end_time+1, time_step):
        draw_topography(time, output_dir, topo_dir, topo_suffix)

