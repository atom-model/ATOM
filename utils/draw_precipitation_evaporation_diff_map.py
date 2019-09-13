#!/usr/bin/python
# -*- coding: utf-8 -*-

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import LinearSegmentedColormap

import os
import numpy as np
from netCDF4 import Dataset
from common_functions import calculate_spherical_mean, convert_atom_to_gmt, add_lon_lat_to_gmt_data

def draw_precipitation_evaporation_diff_map(time=0, data_dir='./', 
        output_dir='./precipitation_evaporation_diff', topo_suffix='smooth'):
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
    tmp_dir = '/tmp/atom/'
    if not os.path.isdir(tmp_dir):
        os.mkdir(tmp_dir)
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)

    np.savetxt(tmp_dir + 'p.xyz',data, fmt='%0.2f')

    os.system(gmt_cmd + ' grdfilter {0}p.xyz -G{0}p.nc -Fm7 -Dp -Vl -fg'.format(tmp_dir))

    fh = Dataset(tmp_dir + 'p.nc', mode='r')
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
    data = evapor_data - data #evaporation minus precipitation

    mean_val = calculate_spherical_mean(data[~np.isnan(data)].flatten(), yy[~np.isnan(data)].flatten())
    print(mean_val)
    topo= all_data[:,2]
    
    plt.figure(figsize=(15, 8))
    m = Basemap(llcrnrlon=-180,llcrnrlat=-90,urcrnrlon=180,urcrnrlat=90,projection='kav7', lon_0=0)

    x, topo = m.shiftdata(x, datain = topo, lon_0=0)

    xi, yi = m(x, y)
    img_data = m.transform_scalar(data, np.arange(-180,180),np.arange(-90,90),361,181)
    cmap = create_cmap()
    cs = m.imshow(img_data,alpha=0.5, vmin=-30*365, vmax=10*365, cmap=cmap)

    m.contour( xi.reshape((361,181)), yi.reshape((361,181)), topo.reshape((361,181)),
                            colors ='k', linewidths= 0.3 )

    m.drawparallels(np.arange(-90., 90., 10.), labels=[1,0,0,0], fontsize=10)
    m.drawmeridians(np.arange(-180., 180., 45.), labels=[0,0,0,1], fontsize=10)

    bounds = np.array([-30,-17,-13,-10,-8,-6,-4 ,-2,-1, -0.2, 0.2, 1, 2 ,4 ,6, 8, 10.0])*365
    norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)
    cb3 = m.colorbar(cs, cmap=cmap,
                                norm=norm,
                                boundaries = bounds,
                                extend='both',
                                extendfrac='auto',
                                ticks=bounds,
                                spacing='uniform',
                                location='bottom',pad="10%", label='Evaporation Minus Precipitation (mm/year)')

    #cbar = m.colorbar(cs, location='bottom', pad="10%", label='Evaporation Minus Precipitation (mm/year)')
    plt.title("{1} at {0}Ma (global mean: {2:.2f})".format(time, 'Evaporation Minus Precipitation', mean_val))
    #plt.annotate('Jul-24-2012', xy=(0.5, 0), xycoords='figure fraction', xytext=(0.5, 0.15), textcoords='figure fraction')
    plt.savefig(output_dir+'/{0}_Ma_{1}.png'.format(time, 'precipitation_evaporation_diff'), bbox_inches='tight')
   
    print(output_dir+'/{0}_Ma_{1}.png has been saved!'.format(time, 'precipitation_evaporation_diff'))
    plt.close()
    fh.close()

def create_cmap():
    cdict = {
            'red':   
                np.array([
                    (0.0,           0.0,        240.0),
                    (13/40.,        240.,       156.0),
                    (17/40.,        156.0,      74.0),
                    (20/40.,        74.0,       16.0),
                    (22/40.,        16.0,       10.0),
                    (24./40.,       10.,        0.0),
                    (29.8/40.,      0.0,        255.0),
                    (1.0,           255.0,      255.0),
            ]),

            'green': 
                np.array([
                    (0.0,           0.0,            0.0),
                    (20/40.,        0.,             20.0),
                    (22/40.,        20.0,           80.0),
                    (24/40.,        80.0,           141.0),
                    (26/40.,        141.0,          155.0),
                    (28./40.,       155.,           145.0),
                    (29/40.,        145.0,          136.0),
                    (29.8/40.,      136.0,          255.0),
                    (31/40.,        255.0,          229.0),
                    (32/40.,        229.0,          204.0),
                    (34/40.,        204.0,          179.0),
                    (36./40.,       179.,           153.0),
                    (38/40.,        153.0,          129.0),
                    (1.0,           129.0,          129.0),

            ]),

            'blue':  
                np.array([
                    (0.0,           0.0,            126.0),
                    (13/40.,        126.,           122.0),
                    (17/40.,        122.0,          124.0),
                    (20/40.,        124.0,          133.0),
                    (22/40.,        133.0,          162.0),
                    (24./40.,       162.,           189.0),
                    (26/40.,        189.0,          163.0),
                    (28/40.,        163.0,          109.0),
                    (29/40.,        109.0,          55.0),
                    (29.8/40.,      55.0,           255.0),
                    (30.2/40.,      255.0,          0.0),
                    (1.,            0.,             0.0),
            ])
    }
    #color_list=[
    #    (0, np.array([255.,129.,0.])/255),
    #    (1, np.array([255.,153.,0.])/255),
    #]
    #return LinearSegmentedColormap.from_list('my_cmap', color_list)
    cdict['red'][:,(1,2)] = cdict['red'][:,(1,2)]/255.
    cdict['green'][:,(1,2)] = cdict['green'][:,(1,2)]/255.
    cdict['blue'][:,(1,2)] = cdict['blue'][:,(1,2)]/255.
    #print(cdict)
    return  LinearSegmentedColormap('my_cmap', segmentdata=cdict)

if __name__ == "__main__":
    for time in range(0,141,10):
        draw_precipitation_evaporation_diff_map(time,'../benchmark/output/')
