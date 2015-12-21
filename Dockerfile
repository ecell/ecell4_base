FROM ubuntu:15.10
RUN apt-get update
RUN apt-get install -y cython git cmake g++ libboost-dev libgsl0-dev libhdf5-dev libboost-regex-dev python python-numpy python-scipy python-pip python-zmq python-matplotlib
RUN pip install jupyter
ADD . /usr/src/ecell4

RUN cd /usr/src/ecell4; export PREFIX=/usr/local; export PYTHONPATH=/usr/local/lib/python2.7/site-packages:$PYTHONPATH; ./install.sh --python2 --hdf5

EXPOSE 8888
CMD LD_LIBRARY_PATH=/usr/local/lib jupyter-notebook --notebook-dir='/usr/src/ecell4/ipynb' --no-browser --ip='*' --port 8888
