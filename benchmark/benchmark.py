#!/usr/bin/env python

from pyatom import Model

model = Model()
model.load_config('benchmark.xml')

# FIXME: there needs to be a trailing slash on output_path; the software should be more tolerant
model.run(t=14, output_path='output-14/')
model.run(t=0, output_path='output-0/')
