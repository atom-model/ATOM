FROM ubuntu:16.04

MAINTAINER Michael Chin

RUN apt-get update -y

RUN apt-get install -y git python-pip gcc

RUN pip install -U setuptools
RUN pip install -U pip
RUN pip install jupyter
RUN pip install netcdf4
RUN pip install colorlover

RUN apt-get install -y wget
# install dependencies for pygplates
RUN apt-get install -y libglew-dev
RUN apt-get install -y python2.7-dev
RUN apt-get install -y libboost-dev libboost-python-dev libboost-thread-dev libboost-program-options-dev libboost-test-dev libboost-system-dev
RUN apt-get install -y libqt4-dev
RUN apt-get install -y libgdal-dev
RUN apt-get install -y libcgal-dev
RUN apt-get install -y libproj-dev
RUN apt-get install -y libqwt-dev
RUN apt-get install -y libxrender-dev libice-dev libsm-dev libfreetype6-dev libfontconfig1-dev

# install wget to enable the pygplates source code to be downloaded from sourceforge
# use wget to get the correct pygplates package from sourceforge
RUN wget http://sourceforge.net/projects/gplates/files/pygplates/beta-revision-12/pygplates-ubuntu-xenial_1.5_1_amd64.deb

# use dpkg to install 
RUN dpkg -i pygplates-ubuntu-xenial_1.5_1_amd64.deb

RUN rm pygplates-ubuntu-xenial_1.5_1_amd64.deb

Env PYTHONPATH ${PYTHONPATH}:/usr/lib:/usr/lib/pygplates/revision12/

# install all the python and ipython notebook requirements
RUN apt-get install -y gcc python-pip
RUN pip install --upgrade pip
RUN pip install numpy scipy matplotlib jupyter pandas sympy nose
RUN pip install ipyparallel pyproj pyshp Pillow
RUN pip install cython

RUN wget https://github.com/matplotlib/basemap/archive/v1.1.0.tar.gz
RUN tar -vxf v1.1.0.tar.gz
RUN cd basemap-1.1.0/ && python setup.py install
RUN rm -rf basemap-1.1.0/ v1.1.0.tar.gz

RUN apt-get install -y gmt python-tk

WORKDIR /home
ARG CACHE_DATE=2019-02-26
RUN git clone https://github.com/atom-model/ATOM.git -b mchin --single-branch
WORKDIR /home/ATOM
RUN make
RUN pip install -e /home/ATOM/python

RUN mkdir /home/ATOM/workspace

# Copy test files to workspace
RUN cp -av /home/ATOM/benchmark/GEOS2115.ipynb /home/ATOM/workspace/
RUN cp -av /home/ATOM/benchmark/GEOS2115_utils.py /home/ATOM/workspace/GEOS2115_utils.py
RUN cp -av /home/ATOM/benchmark/create_maps.ipynb /home/ATOM/workspace/
RUN cp -av /home/ATOM/benchmark/config_atm.xml /home/ATOM/workspace/config_atm.xml
RUN cp -av /home/ATOM/benchmark/config_hyd.xml /home/ATOM/workspace/config_hyd.xml

RUN apt-get update; apt-get install -y curl \
    && curl https://www.earthbyte.org/webdav/ftp/incoming/mchin/Paleotopography_Grids/topo_grids.tar.gz \
    | tar -xvC /home/ATOM/data 

# setup space for working in
VOLUME /home/ATOM/workspace/output

RUN chmod a+x /home/ATOM/docker/notebook.sh

#ENTRYPOINT ["/workspace/notebook.sh"]
CMD ["/home/ATOM/docker/notebook.sh"] 

# Clean up APT when done.
RUN apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

WORKDIR /home/ATOM/workspace

EXPOSE 8888
