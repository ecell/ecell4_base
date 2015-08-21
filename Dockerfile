FROM ubuntu:14.04
RUN apt-get update
RUN apt-get install -y git cmake g++ libboost-dev libgsl0-dev libhdf5-serial-dev libboost-regex-dev python python-numpy python-scipy python-pip python-zmq python-matplotlib
RUN pip install cython "ipython[notebook]"
RUN cd /; git clone git://github.com/ecell/ecell4
RUN cd /ecell4; CPLUS_INCLUDE_PATH=/usr/include cmake .; make; make install
RUN cd /ecell4/python; python setup.py build_ext install

EXPOSE 8888
CMD LD_LIBRARY_PATH=/usr/local/lib ipython notebook --notebook-dir='/ecell4/ipynb' --no-browser --ip='*' --port 8888
