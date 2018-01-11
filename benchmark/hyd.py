#!/usr/bin/env python

import os
import numpy as np
from pyatom import Model

model = Model()
times=range(0,150,10)

for t in range(len(times)):
    time = times[t]
    #print time
    model.run_hyd( time, './output/', './config_hyd.xml' )
