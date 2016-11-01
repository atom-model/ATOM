#!/usr/bin/env python

from pyatom import Atmosphere

model = Atmosphere()
model.load_config('example_atmosphere.xml')
model.run()
