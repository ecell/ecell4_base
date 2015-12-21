FROM ubuntu:15.10
RUN apt-get update
RUN apt-get install -y cython git cmake g++ libboost-dev libgsl0-dev libhdf5-serial-dev libboost-regex-dev python python-numpy python-scipy python-pip python-zmq python-matplotlib
RUN pip install jupyter
ADD . /usr/src/ecell4
RUN cd /usr/src/ecell4; ./install.sh py2

EXPOSE 8888
CMD LD_LIBRARY_PATH=/usr/local/lib jupyter-notebook --notebook-dir='/usr/src/ecell4/ipynb' --no-browser --ip='*' --port 8888
