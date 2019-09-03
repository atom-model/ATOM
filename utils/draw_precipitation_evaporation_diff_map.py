#!/usr/bin/python
# -*- coding: utf-8 -*-

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

import os
import numpy as np
from netCDF4 import Dataset

def calculate_spherical_mean(values, latitudes):
    lat_r=(np.absolute(latitudes))/90.0*np.pi/2
    lat_c= np.cos(np.absolute(lat_r))
    count=0
    for n,l in zip(values, lat_c):
        count+=n*l
    return count/np.sum(lat_c)


def convert_atom_to_gmt(data):
    new_data = np.zeros((181, 361))
    atom_data = data.flatten()
    for i in range(181):
        for j in range(361):
            if j>179:
                j1=j-180
            else:
                j1=j+180
            new_data[i][j] = atom_data[j1*181+i]
    return new_data

def add_lon_lat_to_gmt_data(data):
    lat = np.linspace(90,-90,181)
    lon = np.linspace(-180,180,361)
    X, Y = np.meshgrid(lon,lat)
    data = data.reshape((181,361))
    return np.stack((X, Y, data),axis=-1)

def draw_precipitation_evaporation_diff_map(time=0, data_dir='./', output_dir='/tmp/atom/', topo_suffix='smooth'):
    gmt_cmd = 'gmt' 
    all_data = np.genfromtxt(data_dir + '/[{0}Ma_{1}.xyz]_PlotData_Atm.xyz'.format(time, topo_suffix), skip_header=1)
    data = all_data[:,8]
    evapor_data = all_data[:,10]
    data = convert_atom_to_gmt(data)
    evapor_data = convert_atom_to_gmt(evapor_data)
    
    data = add_lon_lat_to_gmt_data(data)
    data=data.flatten()
    data=data.reshape((len(data)/3,3))
    #print(data)
    np.savetxt('/tmp/atom/p.xyz',data, fmt='%0.2f')

    os.system(gmt_cmd + ' grdfilter /tmp/atom/p.xyz -G/tmp/atom/p.nc -Fm7 -Dp -Vl -fg')

    fh = Dataset('/tmp/atom/p.nc', mode='r')
    #x = fh.variables['lon'][:]
    #y = fh.variables['lat'][:]
    x = all_data[:,0]
    y = all_data[:,1]
    yy = convert_atom_to_gmt(y)
    #x,y=np.meshgrid(x,y)
    data = fh.variables['z'][:]
    data = data*365
    evapor_data = evapor_data*365
    evapor_data = np.flipud(evapor_data)
    #evapor_data[evapor_data==0]=np.nan
    data = data - evapor_data

    mean_val = calculate_spherical_mean(data[~np.isnan(data)].flatten(), yy[~np.isnan(data)].flatten())
    print(mean_val)
    topo= all_data[:,2]
    
    plt.figure(figsize=(15, 8))
    m = Basemap(llcrnrlon=-180,llcrnrlat=-90,urcrnrlon=180,urcrnrlat=90,projection='kav7', lon_0=0)

    x, topo = m.shiftdata(x, datain = topo, lon_0=0)

    xi, yi = m(x, y)
    img_data = m.transform_scalar(data, np.arange(-180,180),np.arange(-90,90),361,181)
    cs = m.imshow(img_data,alpha=0.5, vmin=-2500, vmax=2500, cmap='jet_r')

    m.contour( xi.reshape((361,181)), yi.reshape((361,181)), topo.reshape((361,181)),
                            colors ='k', linewidths= 0.3 )

    m.drawparallels(np.arange(-90., 90., 10.), labels=[1,0,0,0], fontsize=10)
    m.drawmeridians(np.arange(-180., 180., 45.), labels=[0,0,0,1], fontsize=10)

    cbar = m.colorbar(cs, location='bottom', pad="10%", label='Difference between Precipitation and Evaporation')
    plt.title("{1} at {0}Ma (global mean: {2:.2f})".format(time, 'Difference between Precipitation and Evaporation', mean_val))
    #plt.annotate('Jul-24-2012', xy=(0.5, 0), xycoords='figure fraction', xytext=(0.5, 0.15), textcoords='figure fraction')
    plt.savefig(output_dir+'/{0}_Ma_{1}.png'.format(time, 'precipitation_evaporation_diff'), bbox_inches='tight')
   
    print(output_dir+'/{0}_Ma_{1}.png has been saved!'.format(time, 'precipitation_evaporation_diff'))
    plt.close()
    fh.close()

if __name__ == "__main__":
    for time in range(0,141,10):
        draw_precipitation_evaporation_ratio_map(time,'../benchmark/output/')
