#!/usr/bin/env python
import os
import numpy as np
from pyatom import Model, Atmosphere, Hydrosphere

model = Model()
times=range(0,150,10)

for t in range(len(times)):
    time = times[t]
    model.run_atm( time, './output/', './config_atm.xml' )
