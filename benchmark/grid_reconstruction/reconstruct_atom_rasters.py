import os
import math
import numpy as np
import points_in_polygons
import points_spatial_tree
from proximity_query import find_closest_geometries_to_points_using_points_spatial_tree
import pygplates

#import inpaint
from netCDF4 import Dataset
import scipy.interpolate as spi



def write_xyz_file(output_filename, output_data):
    with open(output_filename, 'w') as output_file:
        for output_line in output_data:
            output_file.write(' '.join(str(item) for item in output_line) + '\n')


def xyzfile_to_spatial_tree_of_points(xyzfile):

    data = np.loadtxt(xyzfile)

    points = [pygplates.PointOnSphere(lat, lon) for lat, lon in zip(data[:,1],data[:,0])]

    return points,data


def reconstruct_raster_stage(static_polygon_features, 
                             rotation_model,
                             time_from,
                             time_to,
                             uniform_recon_points, 
                             spatial_tree_of_uniform_recon_points):
    
    print 'Reconstruct static polygons...'
    
    # Reconstruct the multipoint feature.
    recon_static_polygon_features = []
    pygplates.reconstruct(static_polygon_features, rotation_model, recon_static_polygon_features, time_to)
    
    # Extract the polygons and plate IDs from the reconstructed static polygons.
    recon_static_polygons = []
    recon_static_polygon_plate_ids = []
    for recon_static_polygon_feature in recon_static_polygon_features:
        recon_plate_id = recon_static_polygon_feature.get_feature().get_reconstruction_plate_id()
        recon_polygon = recon_static_polygon_feature.get_reconstructed_geometry()
        
        recon_static_polygon_plate_ids.append(recon_plate_id)
        recon_static_polygons.append(recon_polygon)
    
    print 'Find static polygons...'
    
    # Find the reconstructed static polygon (plate IDs) containing the uniform (reconstructed) points.
    #
    # The order (and length) of 'recon_point_plate_ids' matches the order (and length) of 'uniform_recon_points'.
    # Points outside all static polygons return a value of None.
    recon_point_plate_ids = points_in_polygons.find_polygons_using_points_spatial_tree(
            uniform_recon_points, spatial_tree_of_uniform_recon_points, recon_static_polygons, recon_static_polygon_plate_ids)
    
    print 'Group by polygons...'
    
    # Group recon points with plate IDs so we can later create one multipoint per plate.
    recon_points_grouped_by_plate_id = {}
    for point_index, point_plate_id in enumerate(recon_point_plate_ids):
        # Reject any points outside all reconstructed static polygons.
        if point_plate_id is None:
            continue
        
        # Add empty list to dict if first time encountering plate ID.
        if point_plate_id not in recon_points_grouped_by_plate_id:
            recon_points_grouped_by_plate_id[point_plate_id] = []
        
        # Add to list of points associated with plate ID.
        recon_point = uniform_recon_points[point_index]
        recon_points_grouped_by_plate_id[point_plate_id].append(recon_point)
    
    print 'Reverse reconstruct points...'
    
    # Reconstructed points.
    recon_point_lons = []
    recon_point_lats = []
    
    # Present day points associated with reconstructed points.
    point_lons = []
    point_lats = []
    
    # Create a multipoint feature for each plate ID and reverse-reconstruct it to get present-day points.
    #
    # Iterate over key/value pairs in dictionary.
    for plate_id, recon_points_in_plate in recon_points_grouped_by_plate_id.iteritems():
        # Reverse reconstructing a multipoint is much faster than individually reverse-reconstructing points.
        multipoint_feature = pygplates.Feature()
        multipoint_feature.set_geometry(pygplates.MultiPointOnSphere(recon_points_in_plate))
        multipoint_feature.set_reconstruction_plate_id(plate_id)
        
        # Reverse reconstruct the multipoint feature.
        pygplates.reverse_reconstruct(multipoint_feature, rotation_model, time_to)

        #Forward reconstruct multipoint to 
        multipoint_at_from_time = []
        pygplates.reconstruct(multipoint_feature,rotation_model,multipoint_at_from_time,time_from)
        
        # Extract reverse-reconstructed geometry.
        multipoint = multipoint_at_from_time[0].get_reconstructed_geometry()
        
        # Collect present day and associated reconstructed points.
        for point_index, point in enumerate(multipoint):
            lat, lon = point.to_lat_lon()
            point_lons.append(lon)
            point_lats.append(lat)
            
            recon_point = recon_points_in_plate[point_index]
            recon_lat, recon_lon = recon_point.to_lat_lon()
            recon_point_lons.append(recon_lon)
            recon_point_lats.append(recon_lat)
    
    print 'Sample present-day grid...'
    
    # Query present-day grid using present-day points.
    #
    # TODO: Note sure what happens in regions where there's no data in grid (need to ignore those points).
    #data = data_grid.ev(point_lons, point_lats)
    #data = [1.0] * len(recon_point_lons)
    #data = sample_grid_using_scipy(point_lons,point_lats,grdfile)
    
    return recon_point_lons,recon_point_lats,point_lons,point_lats


def run_grid_pip(recon_time, points, polygons, rotation_model):

    reconstructed_polygons = []
    pygplates.reconstruct(polygons, rotation_model, reconstructed_polygons, recon_time)
    rpolygons = []
    for polygon in reconstructed_polygons:
        if polygon.get_reconstructed_geometry():
            rpolygons.append(polygon.get_reconstructed_geometry())
    polygons_containing_points = points_in_polygons.find_polygons(points, rpolygons)
    lat = []
    lon = []
    zval = []
    for pcp,point in zip(polygons_containing_points,points):
        lat.append(point.get_latitude())
        lon.append(point.get_longitude())
        if pcp is not None:
            zval.append(1)
        else:
            zval.append(0)
    return zval


def run_grid_pnp(recon_time, points, spatial_tree_of_uniform_recon_points, polygons, rotation_model, distance_threshold_radians=2):

    reconstructed_polygons = []
    pygplates.reconstruct(polygons, rotation_model, reconstructed_polygons, recon_time)
    rpolygons = []
    for polygon in reconstructed_polygons:
        if polygon.get_reconstructed_geometry():
            rpolygons.append(polygon.get_reconstructed_geometry())
    res = find_closest_geometries_to_points_using_points_spatial_tree(points,
                                                                     spatial_tree_of_uniform_recon_points,
                                                                     rpolygons,
                                                                     distance_threshold_radians = distance_threshold_radians,
                                                                     geometries_are_solid = True)

    pnp_test = []
    for index in res:
        if index is not None:
            pnp_test.append(1)
        else:
            pnp_test.append(0)

    return pnp_test


###########################################################

