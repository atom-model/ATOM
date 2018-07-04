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

WORKDIR /build
ARG CACHE_DATE=2016-01-01
RUN git clone https://github.com/atom-model/ATOM.git -b Roger --single-branch
WORKDIR /build/ATOM
RUN make
RUN pip install -e /build/ATOM/python

RUN mkdir /workspace

# Copy test files to workspace
RUN cp -av /build/ATOM/examples/* /workspace/

# setup space for working in
VOLUME /workspace/volume

RUN chmod a+x /build/ATOM/docker/notebook.sh

#ADD notebook.sh /workspace/
#RUN chmod a+x /workspace/notebook.sh
#ENTRYPOINT ["/workspace/notebook.sh"]
CMD ["/build/ATOM/docker/notebook.sh"] 

# Clean up APT when done.
RUN apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

WORKDIR /build/ATOM/benchmark

EXPOSE 8888
