import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np
import math, os

for t in range(0,231):
    os.system('gmt grd2xyz paleotopobathy_smooth_{0}.00Ma.nc | gmt surface -Rd -I1 -Gpaleotopobathy_smooth_{0}.00Ma_new.nc'.format(t))
    os.system('gmt grd2xyz paleotopobathy_smooth_{0}.00Ma_new.nc > {0}Ma_Simon.xyz'.format(t))
