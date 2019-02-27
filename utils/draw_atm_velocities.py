#!/usr/bin/python
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

def down_sample(a):
    tmp=[]   
    for t in a:
        tmp.append(t[3::4])
    return tmp[3::4]

def draw_velocity(time, output_dir, data_dir):
    v = np.fromfile('./output/bin_data/v_{}_n_0.bin'.format(time),'<f8')
    w = np.fromfile('./output/bin_data/w_{}_n_0.bin'.format(time),'<f8')
    h = np.fromfile('./output/bin_data/h_{}_n_0.bin'.format(time),'<f8') 

    vm = np.sqrt(v**2+w**2)
    vm[vm==0] = 1
    v = v/vm
    w = w/vm

    x = np.linspace(-180, 180, 361)
    y = np.linspace(-90, 90, 181)
    xv, yv = np.meshgrid(x, y)
    #print xv
    #print yv

    figure = plt.figure(figsize=(15, 8))
    m = Basemap(llcrnrlon=-180,llcrnrlat=-90,urcrnrlon=180,urcrnrlat=90,projection='cyl', lon_0=0)

    xi, yi = m(xv.flatten(), yv.flatten())

    xi = xi.reshape((181,361)) 
    yi = yi.reshape((181,361)) 
    vm = vm.reshape((181,361)) 
    w = w.reshape((181,361)) 
    v = v.reshape((181,361)) 
    h = h.reshape((181,361)) 

    
    wn = down_sample(w) 
    vn = down_sample(-v)
    wn = wn / np.sqrt(np.square(wn) + np.square(vn))
    vn = vn / np.sqrt(np.square(wn) + np.square(vn))
    cs = m.quiver(down_sample(xi), down_sample(yi), wn, vn, down_sample(vm), width=0.001,
         headlength=7, headwidth=5, pivot='tail', clim=[0, 1.2], cmap='jet')


    m.contour( xi, yi, h, colors ='k', linewidths= 0.3 )
    m.drawparallels(np.arange(-90., 90., 10.), labels=[1,0,0,0], fontsize=10)
    m.drawmeridians(np.arange(-180., 180., 45.), labels=[0,0,0,1], fontsize=10)

    cbar = m.colorbar(cs, location='bottom', pad="10%", label='Velocity (m/s)')
    plt.title("Atmospheric Velocity at {0}Ma".format(time))
    plt.savefig(output_dir + '/{}Ma_velocity.png'.format(time), bbox_inches='tight')
    plt.close()

if  __name__ == "__main__":
    draw_velocity(0, '../benchmark/output/atm_maps/velocity/', '../benchmark/output/')
    draw_velocity(5, '../benchmark/output/atm_maps/velocity/', '../benchmark/output/')
