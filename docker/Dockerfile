# Pull base image.
FROM phusion:baseimage

MAINTAINER Ian Howson

# Use baseimage-docker's init system.
CMD ["/sbin/my_init"]

#RUN apt-get update -y

#RUN apt-get install -y git python-pip python-dev libzmq3 libzmq3-dev pkg-config libfreetype6-dev libpng3 libopenmpi-dev openmpi-bin libhdf5-dev liblapack-dev llvm-3.6 libedit-dev gfortran

#RUN pip install -U setuptools
#RUN pip install -U pip  # fixes AssertionError in Ubuntu pip
#RUN pip install enum34
#RUN LLVM_CONFIG=llvm-config-3.6 pip install llvmlite==0.8.0
#RUN pip install jupyter markupsafe zmq singledispatch backports_abc certifi jsonschema ipyparallel path.py matplotlib mpi4py==1.3.1 git+https://github.com/badlands-model/triangle pandas plotly
#RUN apt-get install -y libnetcdf-dev python-mpltoolkits.basemap
#RUN pip install Cython==0.20
#RUN pip install h5py
#RUN pip install scipy
#RUN pip install numpy
#RUN pip install numba==0.23.1 ez_setup
#RUN pip install gFlex
#RUN pip install netcdf4
#RUN pip install colorlover
#RUN pip install pyevtk

WORKDIR /build
# FIXME need a public git repo
RUN git clone https://github.com/badlands-model/pyBadlands.git
WORKDIR /build/pyBadlands/pyBadlands/libUtils
RUN make
RUN pip install -e /build/pyBadlands

WORKDIR /build
RUN git clone https://github.com/badlands-model/pyBadlands-Companion.git
RUN pip install -e /build/pyBadlands-Companion

ENV TINI_VERSION v0.8.4
ADD https://github.com/krallin/tini/releases/download/${TINI_VERSION}/tini /usr/local/bin/tini
RUN chmod +x /usr/local/bin/tini

RUN mkdir /workspace && \
    mkdir /workspace/volume && \
        mkdir /workspace/companion

# Copy test files to workspace
RUN cp -av /build/pyBadlands/Examples/* /workspace/
RUN cp -av /build/pyBadlands-Companion/notebooks/* /workspace/companion/

COPY run.sh /build
RUN chmod +x /build/run.sh

# setup space for working in
VOLUME /workspace/volume

# launch notebook
WORKDIR /workspace
EXPOSE 8888
ENTRYPOINT ["/usr/local/bin/tini", "--"]

#ENV LD_LIBRARY_PATH=/workspace/volume/pyBadlands/pyBadlands/libUtils:/build/pyBadlands/pyBadlands/libUtils
#CMD /build/run.sh


RUN mkdir /etc/service/memcached
ADD memcached.sh /etc/service/memcached/run

# Clean up APT when done.
RUN apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# see for more information https://github.com/phusion/baseimage-docker#using
