#!/usr/bin/python
# -*- coding: utf-8 -*-

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import mpl_toolkits
mpl_toolkits.__path__.append('/usr/local/lib/python2.7/dist-packages/mpl_toolkits/')
from mpl_toolkits.basemap import Basemap

import os
import numpy as np
from netCDF4 import Dataset

from common_functions import calculate_spherical_mean

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

def draw_precipitation_map(time=0, data_dir='./', output_dir='/tmp/atom/', topo_suffix='smooth'):
    gmt_cmd = 'gmt' 
    all_data = np.genfromtxt(data_dir + '/[{0}Ma_{1}.xyz]_PlotData_Atm.xyz'.format(time, topo_suffix), skip_header=1)
    data = all_data[:,8]
    data = convert_atom_to_gmt(data)
    data = add_lon_lat_to_gmt_data(data)
    data=data.flatten()
    data=data.reshape((len(data)/3,3))
    #print(data)
    tmp_dir = '/tmp/atom/'
    if not os.path.isdir(tmp_dir):
        os.mkdir(tmp_dir)
    np.savetxt(tmp_dir+'p.xyz',data, fmt='%0.2f')

    os.system(gmt_cmd + ' grdfilter {0}p.xyz -G{0}p.nc -Fm7 -Dp -Vl -fg'.format(tmp_dir))

    fh = Dataset(tmp_dir+'p.nc', mode='r')
    #x = fh.variables['lon'][:]
    yy = fh.variables['lat'][:]
    x = all_data[:,0]
    y = all_data[:,1]
    #x,y=np.meshgrid(x,y)
    data=fh.variables['z'][:]
    data=data*365

    yyy = np.array([[n]*361 for n in yy]).flatten()
    data_1 = data.flatten()
    mean_val = calculate_spherical_mean(data_1[~np.isnan(data_1)], yyy[~np.isnan(data_1)])

    topo= all_data[:,2]
    
    plt.figure(figsize=(15, 8))
    m = Basemap(llcrnrlon=-180,llcrnrlat=-90,urcrnrlon=180,urcrnrlat=90,projection='kav7', lon_0=0)

    x, topo = m.shiftdata(x, datain = topo, lon_0=0)

    xi, yi = m(x, y)
    img_data = m.transform_scalar(data, np.arange(-180,180),np.arange(-90,90),361,181)
    cs = m.imshow(img_data,alpha=0.5, vmin=0, vmax=2, cmap='jet')

    m.contour( xi.reshape((361,181)), yi.reshape((361,181)), topo.reshape((361,181)),
                            colors ='k', linewidths= 0.3 )

    m.drawparallels(np.arange(-90., 90., 10.), labels=[1,0,0,0], fontsize=10)
    m.drawmeridians(np.arange(-180., 180., 45.), labels=[0,0,0,1], fontsize=10)

    cbar = m.colorbar(cs, location='bottom', pad="10%", label='Precipitation (mm/yr)')
    plt.title("{1} at {0}Ma (global spherical mean: {2:.2f})".format(time, 'Precipitation', mean_val))
    plt.savefig(output_dir+'/{0}_Ma_{1}.png'.format(time, 'precipitation'), bbox_inches='tight')
    #print(output_dir+'/precipitation'+'/{0}_Ma_{1}.png has been saved!'.format(time, 'precipitation'))
    plt.close()
    fh.close()

def draw_nasa_precipitation_and_diff_map(time=0, data_dir='./', output_dir='/tmp/atom/', topo_suffix='smooth'):
    data = np.genfromtxt('../data/SurfacePrecipitation_NASA.xyz')
    topo_data = np.genfromtxt('../data/topo_grids/0Ma_smooth.xyz')

    x = data[:,0]
    y = data[:,1]
    z = data[:,2] * 365
    xt = topo_data[:,0]
    yt = topo_data[:,1]
    topo = topo_data[:,2]
    topo = np.where(topo>0,1, 0)

    fig = plt.figure(figsize=(15, 8))
    m = Basemap(llcrnrlon=-180,llcrnrlat=-90,urcrnrlon=180,urcrnrlat=90,projection='kav7', lon_0=0)

    x, z = m.shiftdata(x, datain = z, lon_0=0)

    xi, yi = m(x, y)
    xit, yit = m(xt, yt)

    zz = np.zeros((181, 361))
    for i in range(181):
        for j in range(361):
            zz[180-i][j] = z[j*181+i]
    img_data = m.transform_scalar(zz, np.arange(-180,180),np.arange(-90,90),361,181)
    cs = m.imshow(img_data,alpha=0.5, clim=[0,4000], cmap='jet')

    m.contour( xit.reshape((181,361)), yit.reshape((181,361)), topo.reshape((181,361)),
                    colors ='k', linewidths= 0.3 )
    m.drawparallels(np.arange(-90., 90., 10.), labels=[1,0,0,0], fontsize=10)
    m.drawmeridians(np.arange(-180., 180., 45.), labels=[0,0,0,1], fontsize=10)

    cbar = m.colorbar(cs, location='bottom', pad="10%", label='(mm/yr)')

    plt.title("NASA Present Day Precipitation")

    plt.savefig(output_dir+'/NASA_present_day_precipitation.png', bbox_inches='tight')
    #plt.show()

    zz=z

    data = np.genfromtxt(data_dir + '/[{0}Ma_{1}.xyz]_PlotData_Atm.xyz'.format(time, topo_suffix), skip_header=1)        
    
    x = data[:,0]
    y = data[:,1]
    z = data[:,8]
    topo = data[:,2]

    z=z*365 #365 days in one year
    z=z-zz
    fig = plt.figure(figsize=(15, 8))
    m = Basemap(llcrnrlon=-180,llcrnrlat=-90,urcrnrlon=180,urcrnrlat=90,projection='kav7', lon_0=0)

    xx = x
    x, topo = m.shiftdata(xx, datain = topo, lon_0=0)
    x, z = m.shiftdata(xx, datain = z, lon_0=0)

    xi, yi = m(x, y)
    v_min = 0
    v_max = 1200     

    zz = np.zeros((181, 361))
    for i in range(181):
        for j in range(361):
            zz[180-i][j] = z[j*181+i]
    img_data = m.transform_scalar(zz, np.arange(-180,180),np.arange(-90,90),361,181)
    cs = m.imshow(img_data,alpha=0.5, vmin=0, vmax=3000, cmap='jet')

    m.contour( xi.reshape((361,181)), yi.reshape((361,181)), topo.reshape((361,181)),
                    colors ='k', linewidths= 0.3 )
    m.drawparallels(np.arange(-90., 90., 10.), labels=[1,0,0,0], fontsize=10)
    m.drawmeridians(np.arange(-180., 180., 45.), labels=[0,0,0,1], fontsize=10)

    cbar = m.colorbar(cs, location='bottom', pad="10%", label='Precipitation (mm/yr)')
    plt.title('ATOM Minus NASA Present Day Precipitation')
    plt.savefig(output_dir+'/ATOM_Minus_NASA_present_day_precipitation.png', bbox_inches='tight')
    #plt.show()

if __name__ == "__main__":
    draw_nasa_precipitation_and_diff_map(0,'../benchmark/output/')
    for time in range(0,141,10):    
        draw_precipitation_map(time,'../benchmark/output/')
