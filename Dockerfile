FROM ubuntu:14.04
RUN apt-get update

# RUN apt-get install -y g++ python-dev libboost-dev libgsl0-dev libhdf5-serial-dev pkg-config cython python-pip git libboost-test-dev python-zmq python-matplotlib libboost-regex-dev
# RUN pip install ipython jinja2 tornado
# RUN cd /; git clone git://github.com/ecell/ecell4
# RUN cd /ecell4; PREFIX=/usr/local PYTHONPATH=/usr/local/lib/python2.7/dist-packages bash /ecell4/install.sh core core_python gillespie gillespie_python ode ode_python lattice lattice_python util_python bd bd_python

RUN apt-get install -y git cmake g++ libboost-dev libgsl0-dev libhdf5-serial-dev libboost-regex-dev python python-numpy cython python-pip python-zmq python-matplotlib
RUN cd /; git clone git://github.com/ecell/ecell4
RUN cd /ecell4; CPLUS_INCLUDE_PATH=/usr/include cmake .; make; make install

RUN cd /ecell4/python; python setup.py install

RUN cd /; pip install ipython jinja2 tornado

EXPOSE 8888
CMD LD_LIBRARY_PATH=/usr/local/lib ipython notebook --notebook-dir='/ecell4/ipynb' --no-browser --ip='*' --port 8888
