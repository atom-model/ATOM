import numpy as np

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

