#!/usr/bin/python
# -*- coding: utf-8 -*-
import os, sys
import pygplates
import numpy as np
from scipy.interpolate import griddata
from scipy.ndimage import gaussian_filter

import matplotlib
matplotlib.use('Agg')
from mpl_toolkits.basemap import Basemap
from shapely import geometry

ROTATION_DIR = os.path.dirname(os.path.realpath(__file__)) + '/data'
DATA_DIR = './output'
BATHYMETRY_SUFFIX = 'Ma_smooth.xyz'

script_dir = os.path.dirname(os.path.realpath(__file__)) 
ATOM_HOME = script_dir + '/../'
TOPO_DIR = ATOM_HOME + '/data/topo_grids/'

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

def convert_gmt_to_atom(data):
    new_data = np.zeros((361, 181))
    for lon in range(0,361):
        for lat in range(90,-91,-1):
            if lon<180:
                l=lon+180
            else:
                l=lon-180
            new_data[lon][90-lat] = data[90-lat][l]
    return new_data

def add_lon_lat_to_gmt_data(data):
    lat = np.linspace(90,-90,181)
    lon = np.linspace(-180,180,361)
    X, Y = np.meshgrid(lon,lat)
    data = data.reshape((181,361))
    return np.stack((X, Y, data),axis=-1)

def add_lon_lat_to_atom_data(data):
    lat = np.linspace(90,-90,181)
    lon = np.linspace(0,360,361)
    Y, X = np.meshgrid(lat,lon)
    return np.stack((X, Y, data),axis=-1)

def interp_grid(a):
    d = np.array(a)
    x=d[:,0]
    y=d[:,1]
    z=d[:,2]
    lat = np.linspace(90,-90,181)
    lon = np.linspace(-180,180,361)
    X, Y = np.meshgrid(lon,lat)
    
    grid_data = griddata((x, y), z, (X, Y), method='cubic', fill_value=0)

    zmax = z.max()
    zmin = z.min()
    grid_data[grid_data>zmax] = zmax
    grid_data[grid_data<zmin] = zmin
    
    return grid_data

def interp_grid_gmt(a):
    data = np.array(a)
    gmt_cmd = 'gmt'
    if not os.path.isdir("/tmp/atom/"):
        os.mkdir('/tmp/atom/')
    with open('/tmp/atom/no_land.xyz', 'w') as of:
        for line in data:
            of.write(' '.join(str(l) for l in line) + '\n')
    
    zmax = data[:,2].max()
    zmin = data[:,2].min()
            
    os.system(gmt_cmd + 
        ' surface /tmp/atom/no_land.xyz -G/tmp/atom/fill_land_gap.nc -Rd -I1. -T0.8 -Ll{0} -Lu{1}'.format(zmin,zmax))        
  
    os.system(gmt_cmd + ' grd2xyz /tmp/atom/fill_land_gap.nc > /tmp/atom/fill_land_gap.xyz')
    
    data = np.genfromtxt('/tmp/atom/fill_land_gap.xyz')   
    
    return data[:,2] 

def gmt_filter(a):
    data = np.array(a)
    gmt_cmd = 'gmt'
    with open('/tmp/atom/result.xyz', 'w') as of:
        for line in data:
            of.write(' '.join(str(l) for l in line) + '\n')
    
    os.system(gmt_cmd + ' grdfilter /tmp/atom/result.xyz -G/tmp/atom/result.nc -Fm7 -Dp -Vl')
    
    os.system(gmt_cmd + ' grd2xyz /tmp/atom/result.nc > /tmp/atom/result_filtered.xyz')
              
    data = np.genfromtxt('/tmp/atom/result_filtered.xyz')   
    
    return data[:,2]

def get_coastline_polygons_from_topography(filename):
    data = np.genfromtxt(filename)   
    m = Basemap(llcrnrlon=-180,llcrnrlat=-90,urcrnrlon=180,urcrnrlat=90,projection='cyl')
    x = data[:,0]
    y = data[:,1]
    topo = data[:,2]
    topo[topo>0]=1
    topo[topo<0]=0
    
    xi, yi = m(x,y)
    xi = xi.reshape((181,361)) 
    yi = yi.reshape((181,361)) 
    topo = topo.reshape((181,361)) 

    cs = m.contourf( xi, yi, topo ,levels=[0,0.9,1.1])
    
    polygons=[]
    for col in cs.collections[1:2]:
        for contour_path in col.get_paths(): 
            for ncp,cp in enumerate(contour_path.to_polygons()):
                x = cp[:,0]
                y = cp[:,1]
                lons, lats = m(x,y,inverse=True)
                new_shape = geometry.Polygon([(i[0], i[1]) for i in zip(lons, lats)])
                if ncp == 0:
                    poly = new_shape
                else:
                    poly = poly.difference(new_shape)
            polygons.append(poly)
 
    polygon_features = []
    for p in polygons:
        x,y = p.exterior.xy
        f = pygplates.Feature()
        f.set_geometry(pygplates.PolygonOnSphere(zip(y,x)))
        polygon_features.append(f)
        
    return polygon_features


def remove_data_nearby_coastline(data, topo, len_to_remove=3):
    data = data.reshape((181,361))
    topo = topo.reshape((181,361))
    for i in range(0,181):
        for j in range(0,361):
            if topo[i][j] and i > 0 and (not topo[i-1][j]): #south coast
                for k in range(1, len_to_remove):
                    if i-k >= 0 and (not topo[i-k][j]):
                        data[i-k][j] = np.nan
                    else:
                        break
            if topo[i][j] and i <180 and (not topo[i+1][j]): #north coast
                for k in range(1, len_to_remove):
                    if i+k <= 180 and (not topo[i+k][j]):
                        data[i+k][j] = np.nan
                    else:
                        break
            if topo[i][j] and j > 0 and (not topo[i][j-1]): #west coast
                for k in range(1, len_to_remove):
                    if j-k >= 0 and (not topo[i][j-k]):
                        data[i][j-k] = np.nan
                    else:
                        break
            if topo[i][j] and j < 360 and (not topo[i][j+1]): #east coast
                for k in range(1, len_to_remove):
                    if j+k < 360 and (not topo[i][j+k]):
                        data[i][j+k] = np.nan
                    else:
                        break
                        

def reconstruct_grid(from_time, input_grid, to_time, output_grid, reconstruction_dir=script_dir):
    print(from_time)
    print(to_time)
    print(input_grid)
    print(output_grid)
    
    data = np.genfromtxt(input_grid)
    data = convert_atom_to_gmt(data[:,2])
    
    if from_time == 0:
        topo = np.genfromtxt(TOPO_DIR+'/0Ma_smooth.xyz')   
        topo = topo[:,2]
        topo[topo>0] = True
        topo[topo<=0] = False
        remove_data_nearby_coastline(data, topo, 5)
        
    data = add_lon_lat_to_gmt_data(data)

    static_polygons = pygplates.FeatureCollection(
        reconstruction_dir+'/data/Muller_etal_AREPS_2016_StaticPolygons.gpmlz' )
        
    rotation_files = [reconstruction_dir + '/data/Rotations/Global_EarthByte_230-0Ma_GK07_AREPS.rot']
    rotation_model = pygplates.RotationModel(rotation_files)

    #use matplotlib contour function to extract polygons from topography data
    coastline_polygons_to_time = get_coastline_polygons_from_topography(TOPO_DIR+'/{}Ma_smooth.xyz'.format(to_time))
    coastline_polygons_from_time = get_coastline_polygons_from_topography(TOPO_DIR+'/{}Ma_smooth.xyz'.format(from_time))
   
    #turn grid data into point feature
    #use subductionZoneDepth and subductionZoneSystemOrder properties to keep grid data and point index
    #the reason of using these two properties is I don't know how to store data in a feature in other way
    #if you know better way to store data in a feature, you may change the code below
    points_in_ocean = []
    points_on_land = []
    data = data.reshape((181*361,3))
    for idx, point in enumerate(data):        
        on_land = False
        for p in coastline_polygons_from_time: 
            if p.get_geometry().is_point_in_polygon((float(point[1]),float(point[0]))):
                f = pygplates.Feature()
                f.set_geometry(pygplates.PointOnSphere(float(point[1]),float(point[0])))
                f.set_double(
                    pygplates.PropertyName.create_gpml('subductionZoneDepth'),
                    point[2])
                f.set_integer(
                    pygplates.PropertyName.create_gpml('subductionZoneSystemOrder'),
                    idx)
                points_on_land.append(f)
                on_land=True
                break
        if not on_land and (not np.isnan(point[2])):
            points_in_ocean.append(point)


    #assign plate id to the points_on_land features
    #we only reconstruct the points on continents
    print('assigning plate ids...')
    partitioned_features, unpartitioned_features = pygplates.partition_into_plates(
            static_polygons,
            rotation_files,
            points_on_land,
            properties_to_copy = [
                pygplates.PartitionProperty.reconstruction_plate_id],
            reconstruction_time = from_time,
            partition_return = pygplates.PartitionReturn.separate_partitioned_and_unpartitioned
            )

    #use equivalent stage rotation to reconstruct the points and save the data 
    new_points_on_land=[]
    print('reconstructing grid data...')
    fs = sorted(partitioned_features, key=lambda x: x.get_reconstruction_plate_id())
    from itertools import groupby
    for key, group in groupby(fs, lambda x: x.get_reconstruction_plate_id()):
        #print key
        fr = rotation_model.get_rotation(to_time, key)
        for f in group:
            ll = (fr * f.get_geometry()).to_lat_lon()
            v = f.get_double(pygplates.PropertyName.create_gpml('subductionZoneDepth'))
            new_points_on_land.append([ll[1], ll[0], v])

    #gmt_filter(points_in_ocean)
    grid_data = interp_grid_gmt(points_in_ocean)#fill the gaps in the grid
    grid_data = grid_data.reshape((181,361))
    
    points_in_ocean_to_time = []
    data = add_lon_lat_to_gmt_data(grid_data)
    data = data.reshape((181*361,3))
    for idx, point in enumerate(data):        
        on_land = False
        for p in coastline_polygons_to_time: 
            if p.get_geometry().is_point_in_polygon((float(point[1]),float(point[0]))):
                on_land=True
                break
        if not on_land:
            points_in_ocean_to_time.append(point)
    
    new_grid_data = points_in_ocean_to_time
    #only keep points which are inside the to_time polygons. 
    #prevent info on continent leaking into oceans
    for row in new_points_on_land: 
        for f in coastline_polygons_to_time:
            if f.get_geometry().is_point_in_polygon((row[1], row[0])):
                new_grid_data.append(row)
                continue
    
    grid_data = interp_grid_gmt(new_grid_data)#fill the gaps in the grid
    
    #grid_data = add_lon_lat_to_gmt_data(grid_data)
    #grid_data = grid_data.reshape((181*361,3))
    #grid_data = gmt_filter(grid_data)
    
    grid_data = grid_data.reshape((181,361))
    
    #grid_data = gaussian_filter(grid_data, sigma=1.5)
    
    
    output_data = convert_gmt_to_atom(grid_data)
    output_data = add_lon_lat_to_atom_data(output_data)
    output_data = output_data.reshape((361*181,3))

    np.savetxt(output_grid,output_data,fmt='%1.2f') 

def reconstruct_velocity_grid(from_time, input_grid, to_time, output_grid, reconstruction_dir=script_dir):
    print(from_time)
    print(to_time)
    print(input_grid)
    print(output_grid)
    
    data = np.genfromtxt(input_grid)
    data = convert_atom_to_gmt(data[:,2])
    
    #if from_time == 0:
    if True:
        topo = np.genfromtxt('../data/topo_grids/{}Ma_smooth.xyz'.format(from_time))   
        topo = topo[:,2]
        topo[topo>0] = True
        topo[topo<=0] = False
        remove_data_nearby_coastline(data, topo, 5)
        
    data = add_lon_lat_to_gmt_data(data)

    #use matplotlib contour function to extract polygons from topography data
    coastline_polygons_to_time = get_coastline_polygons_from_topography(TOPO_DIR+'/{}Ma_smooth.xyz'.format(to_time))
    coastline_polygons_from_time = get_coastline_polygons_from_topography(TOPO_DIR+'/{}Ma_smooth.xyz'.format(from_time))
   
    #get points in ocean
    points_in_ocean = []
    data = data.reshape((181*361,3))
    for idx, point in enumerate(data):        
        on_land = False
        for p in coastline_polygons_from_time: 
            if p.get_geometry().is_point_in_polygon((float(point[1]),float(point[0]))):
                on_land=True
                break
        if not on_land and (not np.isnan(point[2])):
            points_in_ocean.append(point)

    #gmt_filter(points_in_ocean)

    grid_data = interp_grid_gmt(points_in_ocean)#fill the gaps in the grid
    grid_data = grid_data.reshape((181,361))

    new_grid_data = []
    data = add_lon_lat_to_gmt_data(grid_data)
    data = data.reshape((181*361,3))
    for idx, point in enumerate(data):        
        for p in coastline_polygons_to_time: 
            if p.get_geometry().is_point_in_polygon((float(point[1]),float(point[0]))):
                point[2]=0
                break
        new_grid_data.append(point[2])
    
    #grid_data = add_lon_lat_to_gmt_data(grid_data)
    #grid_data = grid_data.reshape((181*361,3))
    #grid_data = gmt_filter(grid_data)
    
    grid_data = np.array(new_grid_data).reshape((181,361))
    
    #grid_data = gaussian_filter(grid_data, sigma=1.5)
    
    output_data = convert_gmt_to_atom(grid_data)
    output_data = add_lon_lat_to_atom_data(output_data)
    output_data = output_data.reshape((361*181,3))

    np.savetxt(output_grid,output_data,fmt='%1.2f')     

def reconstruct_temperature(time_0, time_1, suffix='Ma_smooth.xyz'):
    st = np.genfromtxt(DATA_DIR + '/[{0}{1}]_PlotData_Atm.xyz'.format(time_0, suffix),skip_header=1)
    data = st[:,[0,1,6]]
    ind = np.lexsort((-data[:,1],data[:,0]))    
    #print(data[ind])
    
    with open(DATA_DIR + '/{0}Ma_Atm_Temperature.xyz'.format(time_0), 'w') as of:
        for l in data[ind]:
            of.write(' '.join(str(item) for item in l) + '\n')

    reconstruct_grid(
        time_0,
        DATA_DIR + '/{0}Ma_Atm_Temperature.xyz'.format(time_0),
        time_1,
        DATA_DIR + '/{0}Ma_Reconstructed_Temperature.xyz'.format(time_1))       


def reconstruct_precipitation(time_0, time_1, suffix='Ma_smooth.xyz'):
    st = np.genfromtxt(DATA_DIR + '/[{0}{1}]_PlotData_Atm.xyz'.format(time_0, suffix),skip_header=1)
    data = st[:,[0,1,8]]
    ind = np.lexsort((-data[:,1],data[:,0]))    
    #print(data[ind])
    
    with open(DATA_DIR + '/{0}Ma_Atm_Precipitation.xyz'.format(time_0), 'w') as of:
        for l in data[ind]:
            of.write(' '.join(str(item) for item in l) + '\n')

    reconstruct_grid(
        time_0,
        DATA_DIR + '/{0}Ma_Atm_Precipitation.xyz'.format(time_0),
        time_1,
        DATA_DIR + '/{0}Ma_Reconstructed_Precipitation.xyz'.format(time_1))  


def reconstruct_salinity(time_0, time_1, suffix='Ma_smooth.xyz'):
    st = np.genfromtxt(DATA_DIR + '/[{0}{1}]_PlotData_Hyd.xyz'.format(time_0, suffix),skip_header=1)
    data = st[:,[0,1,7]]
    ind = np.lexsort((-data[:,1],data[:,0]))
    #print(data[ind])

    with open(DATA_DIR + '/{0}Ma_Hyd_Salinity.xyz'.format(time_0), 'w') as of:
        for l in data[ind]:
            of.write(' '.join(str(item) for item in l) + '\n')

    reconstruct_grid(
        time_0,
        DATA_DIR + '/{0}Ma_Hyd_Salinity.xyz'.format(time_0),
        time_1,
        DATA_DIR + '/{0}Ma_Reconstructed_Salinity.xyz'.format(time_1))


def reconstruct_wind_v(time_0, time_1, suffix='Ma_smooth.xyz'):
    st = np.genfromtxt(DATA_DIR + '/[{0}{1}]_PlotData_Atm.xyz'.format(time_0, suffix),skip_header=1)
    data = st[:,[0,1,3]]
    ind = np.lexsort((-data[:,1],data[:,0]))
    #print(data[ind])

    with open(DATA_DIR + '/{0}Ma_Atm_v.xyz'.format(time_0), 'w') as of:
        for l in data[ind]:
            of.write(' '.join(str(item) for item in l) + '\n')

    reconstruct_velocity_grid(
        time_0,
        DATA_DIR + '/{0}Ma_Atm_v.xyz'.format(time_0),
        time_1,
        DATA_DIR + '/{0}Ma_Reconstructed_wind_v.xyz'.format(time_1))

def reconstruct_wind_w(time_0, time_1, suffix='Ma_smooth.xyz'):
    st = np.genfromtxt(DATA_DIR + '/[{0}{1}]_PlotData_Atm.xyz'.format(time_0, suffix),skip_header=1)
    data = st[:,[0,1,4]]
    ind = np.lexsort((-data[:,1],data[:,0]))
    #print(data[ind])

    with open(DATA_DIR + '/{0}Ma_Atm_w.xyz'.format(time_0), 'w') as of:
        for l in data[ind]:
            of.write(' '.join(str(item) for item in l) + '\n')

    reconstruct_velocity_grid(
        time_0,
        DATA_DIR + '/{0}Ma_Atm_w.xyz'.format(time_0),
        time_1,
        DATA_DIR + '/{0}Ma_Reconstructed_wind_w.xyz'.format(time_1))

def test(filename):
    from_time = 0
    from_file = filename
    for t in range(5,100,5):
        to_file = '{}Ma.xyz'.format(t)
        reconstruct_grid(from_time, from_file, t, to_file)
        from_time = t
        from_file = to_file


def main():
    try:
        if len(sys.argv) == 2:
            test(sys.argv[1])
            return
        global DATA_DIR
        global BATHYMETRY_SUFFIX 
        time_0 = int(sys.argv[1])
        time_1 = int(sys.argv[2])
        DATA_DIR = sys.argv[3]
        BATHYMETRY_SUFFIX = sys.argv[4]
        atm_or_hyd = sys.argv[5]
        #print(time_0)
        #print(time_1)
        if atm_or_hyd == 'atm':
            reconstruct_temperature(time_0, time_1, BATHYMETRY_SUFFIX)
            reconstruct_precipitation(time_0, time_1, BATHYMETRY_SUFFIX)
            reconstruct_wind_v(time_0, time_1, BATHYMETRY_SUFFIX)
            reconstruct_wind_w(time_0, time_1, BATHYMETRY_SUFFIX)
        else:
            reconstruct_salinity(time_0, time_1, BATHYMETRY_SUFFIX)
    except:
        print("Usage: python reconstruct_atom_data.py 0 10 ./output Ma_smooth.xyz atm/hyd") 
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()


