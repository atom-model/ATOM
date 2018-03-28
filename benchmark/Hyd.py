import os
import numpy as np
from pyatom import Model, Atmosphere, Hydrosphere

model = Model()
times=range(0,55,5)

for t in range(len(times)):
    time = times[t]
    model.run_hyd( time, './output/', './config_hyd.xml' )

