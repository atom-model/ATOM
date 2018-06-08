import sys, random
sys.path.append('../reconstruction')
from reconstruct_atom_data import *

from pyatom import Model, Atmosphere, Hydrosphere

model = Model()
times=range(0,55,5)

for t in range(len(times)):
    time = times[t]
    model.run_atm( time, './output/', './config_atm.xml' )
    model.run_hyd( time, './output/', './config_hyd.xml' )
    if t<len(times)-1:
        reconstruct_temperature(time,times[t+1]) 
        reconstruct_precipitation(time,times[t+1])
        reconstruct_salinity(time,times[t+1])
