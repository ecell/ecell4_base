FROM ubuntu:14.04

# E-Cell4 wheel
RUN apt-get update; apt-get install -y python python-dev python3 python3-dev cmake gcc g++ libboost-dev libgsl0-dev libhdf5-dev wget; wget https://bootstrap.pypa.io/get-pip.py; python2 get-pip.py; python3 get-pip.py; pip2 install cython; pip3 install cython
ADD . /usr/src/ecell4
RUN cd /usr/src/ecell4; cmake .; make BesselTables; cd python; python2 setup.py build_ext; python2 setup.py bdist_wheel; python3 setup.py build_ext; python3 setup.py bdist_wheel; pip2 install dist/ecell-4.0.0-cp27-cp27mu-linux_x86_64.whl; pip3 install dist/ecell-4.0.0-cp34-cp34m-linux_x86_64.whl

# matplotlib and jupyter
#RUN apt-get install -y libfreetype6-dev libpng-dev pkg-config python-numpy python3-numpy pandoc
#RUN pip2 install matplotlib jupyter; pip3 install matplotlib jupyter; python2 -m ipykernel install

# ffmpeg and avconv
#RUN apt-get install -y software-properties-common
#RUN add-apt-repository ppa:mc3man/trusty-media -y; apt-get update; apt-get install -y ffmpeg
#RUN apt-get install -y libav-tools

# biomodels sbml parser
#RUN apt-get install -y python-lxml; pip install beautifulsoup4 requests
#RUN pip install python-libsbml
# install with install.sh
#RUN cd /usr/src/ecell4; export PREFIX=/usr/local; export PYTHONPATH=/usr/local/lib/python2.7/site-packages:$PYTHONPATH; ./install.sh --python2 --hdf5

#EXPOSE 8888
#CMD LD_LIBRARY_PATH=/usr/local/lib jupyter-notebook --notebook-dir='/usr/src/ecell4/ipynb' --no-browser --ip='*' --port 8888
#CMD jupyter-notebook --notebook-dir='/usr/src/ecell4/ipynb' --no-browser --ip='*' --port 8888
