#!/usr/bin/env python

from pyatom import Atmosphere

# Create an Atmosphere model with the default configuration
model = Atmosphere()

# If you want a custom configuration, you can load it like so
# model.load_config('example_atmosphere.xml')

# You can also modify model parameters directly through Python
# print model.verbose

# Run the model to a specified time period (in this case, 20Ma)
# model.run(t=20)

# Run the model to completion (specified in the config file)
model.run()
