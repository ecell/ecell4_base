# E-Cell System version 4 

## What is E-Cell System?

E-Cell System, a software platform for modeling, simulation and analysis of complex, heterogeneous and multi-scale systems like the cell.

## Installing E-Cell (Windows)

### Requirements

- Python 2.7.9(**32bit**) https://www.python.org/ftp/python/2.7.9/python-2.7.9.msi
- HDF5-1.8.14 Pre-built Binary(**32-bit**) http://www.hdfgroup.org/ftp/HDF5/current/bin/windows/extra/hdf5-1.8.14-win32-vs2008-shared.zip

Please add `C:\Python27\Scripts` and `C:\Program Files (x86)\HDF_Group\HDF5\1.8.14\bin` to your **PATH** enviromental variable.

And run following command with command prompt.
```
:: Please download ecell4-4.0.0_beta1-cp27-none-win32.whl from release page
pip install ecell4-4.0.0_beta1-cp27-none-win32.whl
```

### IPython notebook
We recommend you run E-Cell4 models from IPython notebook.
Below is IPython notebook(and matplotlib) installation for Windows.

```
pip install python-dateutil
pip install pyparsing
pip install "ipython[notebook]"
```

next, download matplotlib and numpy from

https://downloads.sourceforge.net/project/matplotlib/matplotlib/matplotlib-1.4.2/windows/matplotlib-1.4.2-cp27-none-win32.whl  
http://sourceforge.net/projects/numpy/files/NumPy/1.9.1/numpy-1.9.1-win32-superpack-python2.7.exe/download

and install these exe files.

## Installing E-Cell (Mac OS X)

### Requirements

- homebrew
- hdf5
- pip

```shell
# here we use homebrew to install hdf5, please install hdf5 to /usr/local/lib
brew install homebrew/science/hdf5 --enable-cxx
brew install wget
wget https://bootstrap.pypa.io/get-pip.py
sudo python get-pip.py
# please download whl file from release page(https://github.com/ecell/ecell4/releases)
pip install ecell4-4.0.0_beta1-cp27-none-macosx_10_10_intel.whl --user
# if you are using Mountain Lion
# pip install ecell4-4.0.0b1-cp27-none-macosx_10_9_intel.whl --user
```

### IPython notebook
We recommend you run E-Cell4 models from IPython notebook.
Below is IPython notebook(and matplotlib) installation for Mac.

```shell
pip install matplotlib --user
pip install "ipython[notebook]" --user
cd ~/Library/Python/2.7/ecell4ipynb
PYTHONPATH=~/Library/Python/2.7/lib/python/site-packages/ ipython notebook
```

now you can see IPython notebooks, please open index.ipynb to see E-Cell4 models.

## Installing E-Cell (Ubuntu 14.04)

```shell
$ sudo apt-get install python-pip libgsl0-dev libhdf5-serial-dev libboost-dev
$ sudo pip install http://dev.e-cell.org/downloads/ecell4/ubuntu-trusty-amd64/latest/ecell4-0.0.0-cp27-none-linux_x86_64.whl
```

## Building and installing E-Cell (Ubuntu 14.04)

```shell
# dependent packages
$ sudo apt-get install libgsl0-dev libboost-dev libboost-test-dev libboost-regex-dev libhdf5-serial-dev
$ sudo apt-get instal python-dev cython

$ wget https://github.com/ecell/ecell4/archive/master.zip   
$ unzip master.zip
$ cd ecell4-master
$ PREFIX=/path/to PYTHONPATH=/path/to/lib/python2.7/site-packages ./install.sh
```

## How to use?

```
$ LD_LIBRARY_PATH=/pat/to/lib PYTHONPATH=/path/to/lib/python2.7/site-packages python
Python 2.7.6 (default, Mar 22 2014, 22:59:56) 
[GCC 4.8.2] on linux2
Type "help", "copyright", "credits" or "license" for more information.
>>> from ecell4.core import *
>>> sp = Species("B.A.C")
>>> print sp.serial()
A.B.C
>>> 
```

## Running E-Cell4 with IPython notebook (by docker)

If you use following docker images, you don't need to do OS dependent installation.  
We have already installed E-Cell4 to docker environment.

### boot2docker (Windows or Mac)

The latest version of boot2docker sets up a host only network adaptor which provides access to the container's ports.

```shell
$ boot2docker ssh
######## in boot2docker
docker@boot2docker:~$ docker pull ecell/ecell4:develop
docker@boot2docker:~$ docker run --rm -i -t -p 8888:8888 ecell/ecell4:develop
```

Then you should be able to access the E-Cell4 IPython notebook server using the IP address reported to you using:

```shell
$ boot2docker ip
```

Typically, it is 192.168.59.103, so please open 192.168.59.103:8888 with your favorite browser.
(But it could get changed by Virtualbox's DHCP implementation.)

### Docker (Linux)

```shell
$ sudo docker pull ecell4/ecell4:develop
$ sudo docker run --rm -i -t -p 8888:8888 ecell/ecell4:develop
```

Open localhost:8888 with your favorite browser.


[![Build Status](https://travis-ci.org/ecell/ecell4.svg?branch=master)](https://travis-ci.org/ecell/ecell4)


## IPython notebooks (tutorials) for E-Cell4

Please see http://nbviewer.ipython.org/github/ecell/ecell4/blob/develop/ipynb/index.ipynb
