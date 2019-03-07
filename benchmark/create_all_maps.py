import sys

def create_all_maps(start_time, end_time, time_step, ATOM_HOME, DATA_DIR):
    sys.path.append(ATOM_HOME + 'reconstruction')
    sys.path.append(ATOM_HOME + 'utils')
    import create_atm_maps, create_hyd_maps
    try:
        topo_dir = ATOM_HOME + 'data/topo_grids/'
        topo_suffix = 'smooth'
        atm_map_output_dir = DATA_DIR + 'atm_maps'
        hyd_map_output_dir = DATA_DIR + 'hyd_maps'

        atm_sub_dirs = ['temperature','v_velocity','w_velocity', 'water_vapour', 'precipitation', 
                    'precipitable_water', 'topography', 'velocity']

        create_atm_maps.create_all_maps(atm_sub_dirs, start_time, end_time, time_step, atm_map_output_dir, 
                DATA_DIR, topo_dir, topo_suffix)

        hyd_sub_dirs = ['temperature','v_velocity','w_velocity', 'salinity', 
                'upwelling', 'downwelling', 'velocity']

        create_hyd_maps.create_all_maps(hyd_sub_dirs, start_time, end_time, time_step, hyd_map_output_dir,
                DATA_DIR, topo_dir, topo_suffix)
    except:
        import traceback
        traceback.print_exc() 

if __name__ == "__main__":
    create_all_maps(int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), sys.argv[4], sys.argv[5])
