FROM ubuntu:14.04
RUN apt-get update
RUN apt-get install -y g++ python-dev libboost-dev libgsl0-dev libhdf5-serial-dev pkg-config cython python-pip git libboost-test-dev python-zmq python-matplotlib
RUN pip install ipython jinja2 tornado
RUN cd /; git clone git://github.com/ecell/ecell4
RUN cd /ecell4; git checkout -t origin/fix_srctree_new; CPLUS_INCLUDE_PATH=/usr/include cmake .; make; make install
EXPOSE 8888
CMD LD_LIBRARY_PATH=/usr/local/lib ipython notebook --notebook-dir='/ecell4/ipynb' --no-browser --ip='*' --port 8888
