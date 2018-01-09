
# coding: utf-8

# In[1]:


import sys
sys.path.append('./grid_reconstruction')

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
        data_dir = '../data',
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
    os.system(gmt_cmd+' surface buffer_tmp.xyz -Ggmt_surface_fill_the_gap.nc -Rd -I%0.8f -T0.1' % 1.)

    # create a seperate gmt grids for reconstructed land (with nan's in oceans) and oceans (using
    # surface to fill gaps everwhere). Then blend together with continents taking precedence
    # BUT: gaps within continent will look funny due to filling from oceans (e.g. within Himalayas)
    os.system(gmt_cmd+' xyz2grd buffer_land_tmp.xyz -Gbuffer_land_tmp.nc -Rd -I%0.8f' % 1.)
    os.system(gmt_cmd+' surface buffer_ocean_tmp.xyz -Gfill_ocean_tmp.nc -Rd -I%0.8f -T0.1' % 1.)
    os.system(gmt_cmd+' grdblend buffer_land_tmp.nc fill_ocean_tmp.nc -Ggmt_grdblend_clobber.nc -Cf -Rd -I%0.8f' % 1.)
    #os.system('/opt/gmt5/bin/gmt grdmath buffer_land_tmp.nc fill_ocean_tmp.nc OR = gmt_grdblend_clobber.nc')
    os.system(gmt_cmd+' grdfilter %s -G%s -Fg%0.2f -D4 -Vl' % (
        'gmt_surface_fill_the_gap.nc','gmt_surface_fill_the_gap_filter.nc',1000))
    os.system(gmt_cmd+' grd2xyz gmt_surface_fill_the_gap_filter.nc > {0}'.format(new_ascii_grid_file) )
    
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
    
    '''for line in new_data:
        if line[2] > 40:
            line[2] = 40
        if line[2] < -40:
            line[2] = -40
    '''
    os.system('rm '+ new_ascii_grid_file)       
    write_xyz_file(new_ascii_grid_file, new_data)
    print "Reconstruction done!"
   


# In[2]:


def reconstruct_temperature(time,times):
    st = np.genfromtxt('./output/[{0}Ma_Golonka.xyz]_PlotData_Atm.xyz'.format(time),skip_header=1)
    data = st[:,[0,1,6]]
    for d in data:
        d[1]=90-d[1]
    ind = np.lexsort((-data[:,1],data[:,0]))    
    print data[ind]
    
    with open('./output/[{0}Ma_Golonka.xyz]_PlotData_Atm_temperature.xyz'.format(time), 'w') as of:
        for l in data[ind]:
            of.write(' '.join(str(item) for item in l) + '\n')

    if t < (len(times)-1):
        if not os.path.isdir('./output'.format(times[t+1])):
            os.mkdir('./output'.format(times[t+1]))
        reconstruct_grid(
            time,
           './output/[{0}Ma_Golonka.xyz]_PlotData_Atm_temperature.xyz'.format(time),
            times[t+1],
            './output/{0}Ma_SurfaceTemperature.xyz'.format(times[t+1]))       


# In[ ]:


def reconstruct_precipitation(time,times):
    st = np.genfromtxt('./output/[{0}Ma_Golonka.xyz]_PlotData_Atm.xyz'.format(time),skip_header=1)
    data = st[:,[0,1,8]]
    for d in data:
        d[1]=90-d[1]
    ind = np.lexsort((-data[:,1],data[:,0]))    
    print data[ind]
    
    with open('./output/[{0}Ma_Golonka.xyz]_PlotData_Atm_precipitation.xyz'.format(time), 'w') as of:
        for l in data[ind]:
            of.write(' '.join(str(item) for item in l) + '\n')

    if t < (len(times)-1):
        if not os.path.isdir('./output'.format(times[t+1])):
            os.mkdir('./output'.format(times[t+1]))
        reconstruct_grid(
            time,
           './output/[{0}Ma_Golonka.xyz]_PlotData_Atm_precipitation.xyz'.format(time),
            times[t+1],
            './output/{0}Ma_SurfacePrecipitation.xyz'.format(times[t+1]))  


# In[ ]:


import os
import numpy as np
from pyatom import Model, Atmosphere, Hydrosphere

model = Model()
times=range(0,150,10)

for t in range(len(times)):
    time = times[t]
    model.run_atm( time, './output/', './config_atm.xml' )
    reconstruct_temperature(time,times)
    reconstruct_precipitation(time,times)


