#!/usr/bin/env python

from pyatom import Atmosphere

# config = AtmosphereConfig()
# config.load_from_file('example_atmosphere.xml')

model = Atmosphere()
model.load_config('example_atmosphere.xml')
model.run()

# model.run()

'''
T.x = 15
print(T.x)
print T.Multiply(4, 5)

T = None


'''
