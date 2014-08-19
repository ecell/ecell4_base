# E-Cell System version 4 

## What is E-Cell System?

E-Cell System, a software platform for modeling, simulation and analysis of complex, heterogeneous and multi-scale systems like the cell.

## Requirements

```shell
$ sudo apt-get install libgsl0-dev libboost-dev libboost-test-dev libboost-regex-dev libhdf5-serial-dev
$ sudo apt-get instal python-dev cython
```
## Ubuntu 14.04 LTS (Trusty Tahr) installation

```shell
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

## Running E-Cell4 with IPython notebook

### boot2docker (Windows or Mac)

The latest version of boot2docker sets up a host only network adaptor which provides access to the container's ports.

```shell
$ boot2docker ssh
######## in boot2docker
docker@boot2docker:~$ docker pull ecell/ecell4:latest
docker@boot2docker:~$ docker run --rm -i -t -p 8888:8888 ecell/ecell4:latest
```

Then you should be able to access the E-Cell4 IPython notebook server using the IP address reported to you using:

```shell
$ boot2docker ip
```

Typically, it is 192.168.59.103, so please open 192.168.59.103:8888 with your favorite browser.
(But it could get changed by Virtualbox's DHCP implementation.)

### Docker (Linux)

```shell
$ sudo docker pull ecell4/ecell4:latest
$ sudo docker run --rm -i -t -p 8888:8888 ecell/ecell4:latest
```

Open localhost:8888 with your favorite browser.


[![Build Status](https://travis-ci.org/ecell/ecell4.svg?branch=master)](https://travis-ci.org/ecell/ecell4)


## IPython notebooks (tutorials) for E-Cell4

Please see [this link to index.ipynb](http://nbviewer.ipython.org/github/ecell/ecell4/blob/master/ipynb/index.ipynb)
