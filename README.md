# ATOM

ATOM (Atmospheric and Ocean Model) is a fast climate model.

## Getting started

The easiest way to get started is with the Docker container.

If you're not familiar with Docker, start by downloading and installing the [Docker Toolbox](https://www.docker.com/products/docker-toolbox). Run Kitematic.

From within Kitematic, search for and run the `atom-model/demo` container.

## Repo contents

* `Makefile`: top-level Makefile which builds everything
* `README.md`: you're reading it!
* `atmosphere`: C++ source code for the Atmosphere model
* `cli`: command line interface source code
* `data`: sample data files with initial conditions
* `docker`: Dockerfile and support file to build the demonstration Docker container
* `examples`: demo files for each of the interfaces
* `hydrosphere`: C++ source code for the Hydrosphere model
* `lib`: common files used by both Atmosphere and Hydrosphere
* `python`: source code for the Python interface
* `tinyxml2`: the [TinyXML-2](http://www.grinninglizard.com/tinyxml2/) XML parser

## Compilation

TODO

## CLI usage

The simplest way to run a model is to configure it through the XML file and then run it, batch style. The model will run to completion and you will have output files to examine once it completes.

To get started with this, look at `examples/config.xml` for the default configuration. You can run it with:

    atm examples/config.xml

Model output will be visible in the `output/` directory.

## Jupyter Notebook usage

The Python module can be installed with:

    pip install -e python

Then, start a Jupyter Notebook server:

    jupyter notebook

From within the Jupyter web interface, you can open `examples/Demo.ipynb`. This includes some basic visualisation of the model output.

## Python script usage

You can write Python scripts to manipulate and run the model. See `examples/sample.py` for an example.

Run:

    python examples/sample.py

You will then have output in the `output/` directory to analyse.

## Configuration

Most configuration is done by modifying an XML file. You can also modify parameters with the Python interface.

### XML

Look at `examples/config.xml`. This file describes all of the parameters. You should make a copy and modify it to suit your purposes.

If you don't include a parameter in your XML file, ATOM will use the default value. The defaults are documented in `examples/config.xml`. You might want to include only the modified parameters in your XML file for clarity.

### Python

Look at `examples/sample.py`. It shows how to modify a parameter through Python.

## Requirements:

* Global bathymetry/topography grid in 1° x 1° spacing -- paleotopography/bathymetry grids (Smith et al. 1994; Golonka et al. 1997) between 140 - 0 Ma are included here, created using agegrid rev.210 (or Earthbyte 2013.2.rot)
* Present day surface temperature: included based on NASA
* Present day precipitation: included from NASA
* Present day salinity: included from NASA

## Compilation

### Ubuntu Linux 16.04 LTS

TODO

* For Ubuntu, we recommend ```apt-get install TODO```.

### macOS Sierra (Homebrew)

TODO

* For OS X, we recommend one of

    brew install Cython
    brew install netcdf

(ADVANCED: Cython can also be installed from pip, which is useful if you're running inside a virtualenv.)

## Authors

Roger Grundmann

Papers
