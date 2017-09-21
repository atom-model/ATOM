# -*- coding: utf-8 -*-
# python reconstruct_atom_grid.py 0 '../data/SurfaceTemperature_NASA.xyz' 10 './temp_1.xyz'
import pygplates
import numpy as np
import matplotlib.pyplot as plt
from reconstruct_atom_rasters import *
from proximity_query import *
import points_spatial_tree
from sphere_tools import sampleOnSphere

def reconstruct_grid(
        time_of_existing_grid, 
        ascii_grid_file, 
        time_for_new_grid, 
        new_ascii_grid_file,
        data_dir = './data',
        buffer_degrees = 15.):

    rotation_model=pygplates.RotationModel(
        data_dir+'/Rotations/Matthews_etal_GPC_2016_410-0Ma_GK07.rot' )
    static_polygon_features = pygplates.FeatureCollection(
        data_dir+'/ContinentalPolygons/Matthews_etal_GPC_2016_ContinentalPolygons.gpmlz' )

    tmpX,tmpY = np.meshgrid(np.arange(-180.,181.,1.),np.arange(-90.,91.,1.)) 
    tmpX = tmpX.flatten()
    tmpY = tmpY.flatten()

    points,data_array = xyzfile_to_spatial_tree_of_points(ascii_grid_file)

    spatial_tree_of_uniform_recon_points = points_spatial_tree.PointsSpatialTree(points)

    ##############
    # Continents
    ##############

    # reconstruct points within continent polygons from t_n+time_step to t_n
    recon_point_lons, recon_point_lats, point_lons, point_lats = \
        reconstruct_raster_stage(
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

    #write_xyz_file(new_ascii_grid_file, zip(tmpX, tmpY, interp_total_temp))
    data = zip(tmpX, tmpY, interp_total_temp)
    new_data = []
    
    for x in range(0,361):
        for y in range(90,-91,-1):
            if x<=180:
                new_data.append(data[361*(y+90)+x+180])
            else:
                new_data.append(data[361*(y+90)+x-180])
            
    write_xyz_file(new_ascii_grid_file, new_data)
    print "Reconstruction done!"
    
#reconstruct_grid(0,'../data/SurfaceTemperature_NASA.xyz',10,'./temp.xyz')

import sys
if __name__ == "__main__":
    args=['','','','','./data',15.]
    for i in range(1,len(sys.argv)):
        args[i-1]=sys.argv[i]
        
    reconstruct_grid(int(args[0]),args[1],int(args[2]),args[3],args[4],args[5],)

