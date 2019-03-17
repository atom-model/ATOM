#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys, random

import pygplates
import numpy as np
import matplotlib.pyplot as plt
from reconstruct_atom_rasters import *
from proximity_query import *
import points_spatial_tree
from sphere_tools import sampleOnSphere

#print(os.path.dirname(os.path.realpath(__file__)))

ROTATION_DIR = os.path.dirname(os.path.realpath(__file__)) + '/data'
BUFFER_DEGREES = 15.
DATA_DIR = './output'
BATHYMETRY_SUFFIX = 'Ma_smooth.xyz'

def reconstruct_grid(
        time_of_existing_grid, 
        ascii_grid_file, 
        time_for_new_grid, 
        new_ascii_grid_file,
        rotation_dir = ROTATION_DIR,
        buffer_degrees = BUFFER_DEGREES):

    rotation_model=pygplates.RotationModel(
        rotation_dir+'/Rotations/Global_EarthByte_230-0Ma_GK07_AREPS.rot' )
    static_polygon_features = pygplates.FeatureCollection(
        rotation_dir+'/ContinentalPolygons/Matthews_etal_GPC_2016_ContinentalPolygons.gpmlz' )

    tmpX,tmpY = np.meshgrid(np.arange(-180.,181.,1.),np.arange(-90.,91.,1.)) 
    tmpX = tmpX.flatten()
    tmpY = tmpY.flatten()

    points,data_array = xyzfile_to_spatial_tree_of_points(ascii_grid_file)

    spatial_tree_of_uniform_recon_points = points_spatial_tree.PointsSpatialTree(points)

    ##############
    # Continents
    ##############

    # reconstruct points within continent polygons from t_n+time_step to t_n
    recon_point_lons, recon_point_lats, point_lons, point_lats =         reconstruct_raster_stage(
            static_polygon_features, 
            rotation_model,
            time_of_existing_grid,
            time_for_new_grid,
            points, 
            spatial_tree_of_uniform_recon_points)


    # interpolate the values from the grid at t_n onto the unreconstructed points (where we 
    # already know the point locations at t_n+time_step)
    d,l = sampleOnSphere(data_array[:,1],
                         data_array[:,0],
                         data_array[:,2],
                         np.hstack(point_lats),
                         np.hstack(point_lons),n=4)

    interp_land_temp = data_array[:,2].ravel()[l]


    ##########
    # Oceans
    ##########

    # find points that are not within, and actually not close to (defined by some buffer distance) 
    # the continent polygons at t_n. This results in a set of points far from continents at t_n,
    # but there may be some overlap with points within continents at t_n+time_step if the time step
    # is large and/or the continents are moving quickly
    pnp_test = run_grid_pnp(time_of_existing_grid,
                            points,
                            spatial_tree_of_uniform_recon_points,
                            static_polygon_features, 
                            rotation_model, 
                            np.radians(buffer_degrees))

    ocean_index = np.where(np.array(pnp_test)==0)


    ##########################
    # Merged Land and Oceans
    ##########################

    merged_lons = np.hstack((np.hstack(recon_point_lons),data_array[ocean_index,0].flatten()))
    merged_lats = np.hstack((np.hstack(recon_point_lats),data_array[ocean_index,1].flatten()))
    merged_temp = np.hstack((interp_land_temp,data_array[ocean_index,2].flatten()))

    merged_lons = ((merged_lons+180) % 360) - 180

    d,l = sampleOnSphere(merged_lats,
                         merged_lons,
                         merged_temp,
                         tmpY,
                         tmpX,n=4)

    interp_total_temp = merged_temp.ravel()[l]

    # Write out an xyz file for each scalar type - each file contains (lon, lat, scalar).
    write_xyz_file('tmp.xyz', zip(tmpX, tmpY, interp_total_temp))
    write_xyz_file('buffer_tmp.xyz', zip(merged_lons, merged_lats, merged_temp))
    write_xyz_file('buffer_land_tmp.xyz', zip(np.hstack(recon_point_lons), 
                                              np.hstack(recon_point_lats), 
                                              interp_land_temp))
    write_xyz_file('buffer_ocean_tmp.xyz', zip(data_array[ocean_index,0].flatten(), 
                                               data_array[ocean_index,1].flatten(), 
                                               data_array[ocean_index,2].flatten()))

    gmt_cmd = 'gmt'
    # create a gmt grid that fills gaps using nearest-neighbour interpolation from python
    # (but gaps are filled by interpolation from neighbouring continents as well as oceans)
    os.system(gmt_cmd+' xyz2grd tmp.xyz -Gnearest_neighbour_fill_the_gap.nc -Rd -I%0.8f' % 1.)

    # create a gmt grid that fills gaps using surface 
    # (but gaps are filled by interpolation from neighbouring continents as well as oceans)
    os.system(gmt_cmd+' surface buffer_tmp.xyz -Ggmt_surface_fill_the_gap.nc -Rd -I%0.8f -T0.8' % 1.)

    # create a seperate gmt grids for reconstructed land (with nan's in oceans) and oceans (using
    # surface to fill gaps everwhere). Then blend together with continents taking precedence
    # BUT: gaps within continent will look funny due to filling from oceans (e.g. within Himalayas)
    os.system(gmt_cmd+' xyz2grd buffer_land_tmp.xyz -Gbuffer_land_tmp.nc -Rd -I%0.8f' % 1.)
    os.system(gmt_cmd+' surface buffer_ocean_tmp.xyz -Gfill_ocean_tmp.nc -Rd -I%0.8f -T0.8' % 1.)
    os.system(gmt_cmd+' grdblend buffer_land_tmp.nc fill_ocean_tmp.nc -Ggmt_grdblend_clobber.nc -Cf -Rd -I%0.8f' % 1.)
    #os.system('/opt/gmt5/bin/gmt grdmath buffer_land_tmp.nc fill_ocean_tmp.nc OR = gmt_grdblend_clobber.nc')
    os.system(gmt_cmd+' grdfilter %s -G%s -Fg%0.2f -D4 -Vl' % (
        'gmt_surface_fill_the_gap.nc','gmt_surface_fill_the_gap_filter.nc',1000))
    os.system(gmt_cmd+' grd2xyz gmt_surface_fill_the_gap_filter.nc > {0}'.format(new_ascii_grid_file) )
    
    #os.system(gmt_cmd+' grdfilter %s -G%s -Fg%0.2f -D4 -Vl' % (
    #    'nearest_neighbour_fill_the_gap.nc','nearest_neighbour_fill_the_gap_filter.nc',1000))
    #os.system(gmt_cmd+' grd2xyz nearest_neighbour_fill_the_gap_filter.nc > {0}'.format(new_ascii_grid_file) )    

    #write_xyz_file(new_ascii_grid_file, zip(tmpX, tmpY, interp_total_temp))
    data = np.loadtxt(new_ascii_grid_file)
    new_data = []
    
    for x in range(0,361):
        for y in range(90,-91,-1):
            if x<=180:
                new_data.append(data[361*(90-y)+x+180])
            else:
                l=list(data[361*(90-y)+x-180])
                l[0]+=360
                new_data.append(l)

    #adjust the reconstructed data
    input_data = np.genfromtxt(ascii_grid_file)
    input_min = np.nanmin(input_data[:,2])
    input_max = min(np.nanmax(input_data[:,2]), 40.0)
    #print input_min, input_max

    for line in new_data:
        if line[2] > input_max:
            line[2] = input_max
        elif line[2] < input_min:
            line[2] = input_min
    
    #print new_data
    #print np.nanmin(np.array(new_data)[:,2]), np.nanmax(np.array(new_data)[:,2]) 
    
    os.system('rm '+ new_ascii_grid_file)       
    write_xyz_file(new_ascii_grid_file, new_data)
    print "Reconstruction done!"
   

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


def main():
    try:
        global DATA_DIR
        global BATHYMETRY_SUFFIX 
        time_0 = int(sys.argv[1])
        time_1 = int(sys.argv[2])
        DATA_DIR = sys.argv[3]
        BATHYMETRY_SUFFIX = sys.argv[4]
        atm_or_hyd = sys.argv[5]
        print(time_0)
        print(time_1)
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


