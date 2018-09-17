import numpy as np
data = np.genfromtxt('../data/SurfaceTemperature_NASA.xyz')
lat=data[:,1]
lat_r=(np.absolute(lat))/90.0*np.pi/2
lat_c= np.cos(np.absolute(lat_r))
#print len(data[:,2])
#print np.mean(data[:,2])
print lat_c
count=0
for n,l in zip(data[:,2],lat_c):
    #print n, l
    count+=n*l
print count/np.sum(lat_c)
