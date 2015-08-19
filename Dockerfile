FROM jupyter/notebook

MAINTAINER Kozo Nishida <knishida@riken.jp>


RUN apt-get install -y cmake libboost-dev libgsl0-dev libhdf5-serial-dev libboost-regex-dev python-numpy python-scipy python-matplotlib
RUN pip install cython
RUN cd /; git clone git://github.com/ecell/ecell4
RUN cd /ecell4; CPLUS_INCLUDE_PATH=/usr/include cmake .; make; make install

RUN cd /ecell4/python; python setup.py build_ext install
