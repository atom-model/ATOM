#!/usr/bin/python
# -*- coding: utf-8 -*-
import os, sys
import pygplates
import numpy as np
from scipy.interpolate import griddata

ROTATION_DIR = os.path.dirname(os.path.realpath(__file__)) + '/data'
DATA_DIR = './output'
BATHYMETRY_SUFFIX = 'Ma_smooth.xyz'

script_dir = os.path.dirname(os.path.realpath(__file__)) 

def convert_atom_to_gmt(data):
    new_data = np.zeros((181, 361))
    atom_data = data.flatten()
    for i in range(181):
        for j in range(361):
            if j>179:
                j1=j-180
            else:
                j1=j+180
            new_data[i][j] = atom_data[j1*181+180-i]
    return new_data

def convert_gmt_to_atom(data):
    new_data = np.zeros((361, 181))
    for lon in range(0,361):
        for lat in range(90,-91,-1):
            if lon<180:
                l=lon+180
            else:
                l=lon-180
            new_data[lon][90-lat] = data[lat+90][l]
    return new_data


def reconstruct_grid(from_time, input_grid, to_time, output_grid, reconstruction_dir=script_dir):
    print(from_time)
    print(to_time)
    print(input_grid)
    print(output_grid)
    data = np.genfromtxt(input_grid)
    data = convert_atom_to_gmt(data[:,2])

    lat = np.linspace(-90,90,181)
    lon = np.linspace(-180,180,361)
    X, Y = np.meshgrid(lon,lat)
    data = np.stack((X, Y, data),axis=-1)

    static_polygon_features = pygplates.FeatureCollection(
            reconstruction_dir+'/data/ContinentalPolygons/Matthews_etal_GPC_2016_ContinentalPolygons.gpmlz' )
    rotation_files = [reconstruction_dir + '/data/Rotations/Global_EarthByte_230-0Ma_GK07_AREPS.rot']
    rotation_model = pygplates.RotationModel(rotation_files)

    #reconstruct static continental polygons
    #the polygons will be used to mask grid data and partition the grid
    static_polygon_features_to_time = []
    static_polygon_features_from_time = []
    pygplates.reconstruct(static_polygon_features, rotation_model, static_polygon_features_to_time, to_time)
    pygplates.reconstruct(static_polygon_features, rotation_model, static_polygon_features_from_time, from_time)
    
    print('reconstructing continental polygons...')
    polygon_features_to_time = []
    polygon_features_from_time = []
    for p in static_polygon_features_to_time:
        f = pygplates.Feature()
        f.set_geometry(pygplates.PolygonOnSphere(p.get_reconstructed_geometry()))
        f.set_reconstruction_plate_id(p.get_feature().get_reconstruction_plate_id())
        polygon_features_to_time.append(f)
    #    #print bound.get_feature().get_reconstruction_plate_id()

    for p in static_polygon_features_from_time:
        f = pygplates.Feature()
        f.set_geometry(pygplates.PolygonOnSphere(p.get_reconstructed_geometry()))
        f.set_reconstruction_plate_id(p.get_feature().get_reconstruction_plate_id())
        polygon_features_from_time.append(f)

    #turn grid data into point feature
    #use subductionZoneDepth and subductionZoneSystemOrder properties to keep grid data and point index
    #the reason of using these two properties is I don't know how to store data in a feature in other way
    #if you know better way to store data in a feature, you may change the code below
    points = []
    data = data.reshape((181*361,3))
    for idx, point in enumerate(data):
        f = pygplates.Feature()
        f.set_geometry(pygplates.PointOnSphere(float(point[1]),float(point[0])))
        f.set_double(
            pygplates.PropertyName.create_gpml('subductionZoneDepth'),
            point[2])
        f.set_integer(
            pygplates.PropertyName.create_gpml('subductionZoneSystemOrder'),
            idx)
        points.append(f)

    #assign plate id to the point features
    #we only reconstruct the points on continents
    print('assigning plate ids...')
    assigned_point_features, unpartitioned_features = pygplates.partition_into_plates(
            polygon_features_from_time,
            rotation_files,
            points,
            properties_to_copy = [
                pygplates.PartitionProperty.reconstruction_plate_id],
            partition_return = pygplates.PartitionReturn.separate_partitioned_and_unpartitioned
            )

    #use equivalent stage rotation to reconstruct the points and save the data 
    new_data=[]
    print('reconstructing grid data...')
    fs = sorted(assigned_point_features,key=lambda x: x.get_reconstruction_plate_id())
    from itertools import groupby
    for key, group in groupby(fs, lambda x: x.get_reconstruction_plate_id()):
        fr = rotation_model.get_rotation(to_time, key, from_time)
        for f in group:
            ll = (fr * f.get_geometry()).to_lat_lon()
            v = f.get_double(pygplates.PropertyName.create_gpml('subductionZoneDepth'))
            new_data.append([ll[1], ll[0], v])
    #print len(new_data)   

    #this step is for removing the data on continents later
    partitioned_features, unpartitioned_features = pygplates.partition_into_plates(
            polygon_features_to_time,
            rotation_files,
            points,
            properties_to_copy = [
                pygplates.PartitionProperty.reconstruction_plate_id],
            partition_return = pygplates.PartitionReturn.separate_partitioned_and_unpartitioned
            )

    #mask out the points inside the polygons at from_time and to_time
    print('merging grid data...')
    masked_points = set()
    for f in assigned_point_features:
        i = f.get_integer(pygplates.PropertyName.create_gpml('subductionZoneSystemOrder'))
        masked_points.add(i)

    tmp = []
    for p in points:
        i = p.get_integer(pygplates.PropertyName.create_gpml('subductionZoneSystemOrder'))
        if i not in masked_points:
            ll = p.get_geometry().to_lat_lon()
            v = p.get_double(pygplates.PropertyName.create_gpml('subductionZoneDepth'))
            tmp.append([ll[1], ll[0], v])  

    d = np.array(tmp)
    x=d[:,0]
    y=d[:,1]
    z=d[:,2]
    lat = np.linspace(-90,90,181)
    lon = np.linspace(-180,180,361)
    X, Y = np.meshgrid(lon,lat)
    #print lat, lon
    #print z.max(), z.min()
    grid_data = griddata((x, y), z, (X, Y), method='nearest',fill_value=0)

    zmax = z.max()
    zmin = z.min()
    grid_data[grid_data>zmax] = zmax
    grid_data[grid_data<zmin] = zmin
    
    masked_points = set()
    for f in partitioned_features:
        i = f.get_integer(pygplates.PropertyName.create_gpml('subductionZoneSystemOrder'))
        masked_points.add(i)

    grid_data = np.stack((X, Y, grid_data),axis=-1)

    grid_data = grid_data.reshape((181*361,3))
    print len(new_data)
    #print new_data
    for idx, row in enumerate(grid_data):
        if idx not in masked_points:
            new_data.append(list(row))
    
    d = np.array(new_data)
    print d.shape
    x=d[:,0]
    y=d[:,1]
    z=d[:,2]
    grid_data = griddata((x, y), z, (X, Y), method='nearest',fill_value=0)

    zmax = z.max()
    zmin = z.min()
    grid_data[grid_data>zmax] = zmax
    grid_data[grid_data<zmin] = zmin

    output_data = convert_gmt_to_atom(grid_data)
    lat = np.linspace(90,-90,181)
    lon = np.linspace(0,360,361)
    X, Y = np.meshgrid(lat,lon)
    
    output_data = np.stack((Y, X, output_data),axis=-1)
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

    reconstruct_grid(
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

    reconstruct_grid(
        time_0,
        DATA_DIR + '/{0}Ma_Atm_w.xyz'.format(time_0),
        time_1,
        DATA_DIR + '/{0}Ma_Reconstructed_wind_w.xyz'.format(time_1))

def test(filename):
    from_time = 0
    from_file = filename
    for t in range(10,60,10):
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


