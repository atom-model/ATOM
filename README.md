# ATOM

ATOM (Atmospheric and Ocean Model) is a fast, low-complexity, non-parallel climate model.

## Run ATOM inside the Docker container

ATOM depends on pygplates, gmt, matplotlib and other libraries. It is time-consuming to install all the dependencies in your computer. The easiest way to get started is with the gplates/atom Docker container which the GPlates development team have prepared for the users. 

First, you need to install Docker and familiarise yourself with basic Docker operations. [Docker download](https://www.docker.com/community-edition) 

Run "docker pull gplates/atom" or use [Kitematic](https://kitematic.com/) to download the container from Docker Hub.

Run "docker run -it --rm -p 18888:8888 gplates/atom" to start the Docker container in interactive mode.

Alternatively, you can run "docker run -d -p 18888:8888 gplates/atom" to start the Docker container in daemon mode.

There will be a Jupyter Notebook available at port 18888.

Run "docker exec -it *container_id* /bin/bash" to attach to the Docker container and get a shell.

There is a volume mounted in the container. Run "docker inspect -f "{{ .Mounts }}" *container_id*".

Run "docker ps" to check your container id.

Alternatively, you can use the "-v" flag in "docker run" command line to mount a volume, for example "-v dir_in_your_host:dir_in_container".

If you run "docker run -it --rm -p 18888:8888 gplates/atom /bin/bash", the Jupyter Notebook will not be started. You need to run /build/ATOM/docker/notebook.sh manually. 

### Run Python script

You can write Python scripts to manipulate and run the model. See `benchmark/run.py` for an example.

Go into bechmark directory and run:

        python run.py

You will then have output in the `output/` directory to analyse.

You can run "utils/create_atm_maps.py" to create and save maps in the directory ./atm_maps.

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

### Ubuntu Linux 14.04 LTS or 16.04 LTS

As root:

    apt-get install cython

As a regular user:

    git clone https://github.com/atom-model/ATOM.git
    cd ATOM
    make

### macOS Sierra (Homebrew)

    brew install pkg-config
    pip install cython

    git clone https://github.com/atom-model/ATOM
    cd ATOM
    make
    pip install -e python/

# ATOM requires many dependencies to run. It is highly recommended to use the Docker container to run ATOM.

## CLI usage

The simplest way to run a model is to configure it through the XML file and then run it, batch style. The model will run to completion and you will have output files to examine once it completes.

To get started with this, look at `cli/config_atm.xml` for the default configuration. You can run it with:

    cd cli
    ./atm config_atm.xml

Model output will be visible in the `output/` directory.

## Jupyter Notebook usage

The Python module can be installed with:

    pip install -e python

Then, start a Jupyter Notebook server:

    jupyter notebook

From within the Jupyter web interface, you can open `benchmark/Demo.ipynb`. This includes some basic visualisation of the model output.

## Python script usage

You can write Python scripts to manipulate and run the model. See `benchmark/run.py` for an example.

Run:
    cd benchmark
    python run.py

You will then have output in the `output/` directory to analyse.

## Configuration

Most configuration is done by modifying an XML file. You can also modify parameters with the Python interface.

### XML

Look at `benchmark/config_atm.xml`. This file describes all of the parameters. You should make a copy and modify it to suit your purposes.

If you don't include a parameter in your XML file, ATOM will use the default value. The defaults are documented in `benchmark/config_atm.xml`. You might want to include only the modified parameters in your XML file for clarity.

### Python

Look at `examples/sample.py`. It shows how to modify a parameter through Python.

## Requirements:

* Global bathymetry/topography grid in 1° x 1° spacing -- paleotopography/bathymetry grids (Smith et al. 1994; Golonka et al. 1997) between 140 - 0 Ma are included here, created using agegrid rev.210 (or Earthbyte 2013.2.rot)
* Present day surface temperature: included based on NASA
* Present day precipitation: included from NASA
* Present day salinity: included from NASA



## Parameters

There are many parameters that you can adjust.

**If you are browsing the source code**, the canonical source for parameter documentation is `param.py`. 

**If you are working with an installed distribution**, the canonical source for parameter documentation is `examples/example.xml`.

Much of the code to implement parameters is autogenerated at compile time. This is performed by `param.py`. It generates some `.inc` files (C++) and `.pyi` (Python) files that implement the relevant interfaces.

## Authors

Roger Grundmann, Michael Chin

Papers

Papers used:

An Ice-Water Saturation Adjustment, Wei-Kuo Tao, Joanne Simpson and Michael McCumber, Notes and Correspondence, Monthly Weather Review, Volume 117, January  1989, p 231

A Description of the Nonhydrostatic Regional Model, Part II, Physical Parameterisation, G. Doms et. al.
