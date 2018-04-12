# Pull base image.
FROM phusion/baseimage

MAINTAINER Ian Howson

# Use baseimage-docker's init system.
CMD ["/sbin/my_init"]

RUN apt-get update -y

RUN apt-get install -y git python-pip cython clang libnetcdf-dev python-mpltoolkits.basemap python-numpy python-scipy python-h5py

RUN pip install -U setuptools
RUN pip install -U pip
RUN pip install jupyter
RUN pip install netcdf4
RUN pip install colorlover

WORKDIR /build
#RUN git clone https://github.com/atom-model/ATOM.git
WORKDIR /build/ATOM
RUN make
RUN pip install -e /build/ATOM/python

RUN mkdir /workspace && \
    mkdir /workspace/volume

# Copy test files to workspace
RUN cp -av /build/ATOM/examples/* /workspace/

# setup space for working in
VOLUME /workspace/volume

RUN mkdir /etc/service/notebook
ADD notebook.sh /etc/service/notebook/run
RUN chmod +x /etc/service/notebook/run

# Clean up APT when done.
RUN apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

EXPOSE 8888