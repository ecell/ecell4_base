E-Cell System version 4 
=======================

[![Build Status](https://travis-ci.org/ecell/ecell4.svg?branch=develop)](https://travis-ci.org/ecell/ecell4)

## What is E-Cell System?

E-Cell System, a software platform for modeling, simulation and analysis of complex, heterogeneous and multi-scale systems like the cell.

- [Installation](#installation)
    - [Windows](#windows-installation)
    - [Mac OS X](#mac-os-x-installation)
    - [Ubuntu vivid](#ubuntu-linux-vivid-vervet-installation)
    - [Ubuntu trusty](#ubuntu-linux-trusty-tahr-installation)

Installation
------------

### Windows installation

#### Requirements

Please use 32bit Python, even if you use 64bit Windows.
We don't support 64bit Python

- Python 2.7.10(**32bit**) https://www.python.org/ftp/python/2.7.10/python-2.7.10.msi
- HDF5-1.8.14 Pre-built Binary(**32-bit**) http://www.hdfgroup.org/ftp/HDF5/current/bin/windows/extra/hdf5-1.8.14-win32-vs2008-shared.zip

Please add `C:\Python27`, `C:\Python27\Scripts` and `C:\Program Files (x86)\HDF_Group\HDF5\1.8.14\bin` to your **PATH** enviromental variable.

And run following command with command prompt.
```
pip install https://github.com/ecell/ecell4/releases/download/4.0.0-beta2/ecell4-4.0.0b2-cp27-none-win32.whl
```

#### IPython notebook
We recommend you run E-Cell4 models from IPython notebook.
Below is IPython notebook(and matplotlib) installation for Windows.

- Install [Visual C++ Compiler for Python 2.7](http://www.microsoft.com/en-us/download/details.aspx?id=44266)
- Install IPython notebook and matplotlib

    ```
    pip install "ipython[notebook]"
    pip install matplotlib
    ```

matplotlib depends on numpy. It takes some time to build numpy, please be patient.

### Mac OS X installation

#### Requirements

- homebrew
- gsl
- hdf5
- pip

```shell
brew install homebrew/science/hdf5
brew install gsl
brew install wget
wget https://bootstrap.pypa.io/get-pip.py
sudo python get-pip.py
pip install https://github.com/ecell/ecell4/releases/download/4.0.0-beta2/ecell4-4.0.0b2-cp27-none-macosx_10_10_intel.whl --user
# For Mountain Lion
# pip install https://github.com/ecell/ecell4/releases/download/4.0.0-beta2/ecell4-4.0.0b2-cp27-none-macosx_10_9_intel.whl --user
```

#### IPython notebook
We recommend you run E-Cell4 models from IPython notebook.
Below is IPython notebook(and matplotlib) installation for Mac.

```shell
pip install matplotlib --user
pip install "ipython[notebook]" --user
cd ~/Library/Python/2.7/ecell4ipynb
PYTHONPATH=~/Library/Python/2.7/lib/python/site-packages/ ipython notebook
```

now you can see IPython notebooks, please open index.ipynb to see E-Cell4 models.


### Ubuntu Linux Vivid Vervet installation
#### Python2 series

```shell
# dependent packages
$ sudo apt-get install cmake libgsl0-dev libboost-regex-dev libhdf5-dev cython

$ wget https://github.com/ecell/ecell4/archive/master.zip   
$ unzip master.zip
$ cd ecell4-master
# By default install.sh tries to install E-Cell4 into /usr/local, in this case you need to use sudo.
# In the following command, we install E-Cell4 into $HOME/ecell4. In this case you do NOT need to use sudo.
$ PREFIX=$HOME/ecell4 ./install.sh py2
```

#### Python3 series

```shell
# dependent packages
$ sudo apt-get install cmake libgsl0-dev libboost-regex-dev libhdf5-dev cython3 python3-pip

$ wget https://github.com/ecell/ecell4/archive/master.zip   
$ unzip master.zip
$ cd ecell4-master
# By default install.sh tries to install E-Cell4 into /usr/local, in this case you need to use sudo.
# In the following command, we install E-Cell4 into $HOME/ecell4. In this case you do NOT need to use sudo.
$ PREFIX=$HOME/ecell4 ./install.sh py3
```

### Ubuntu Linux Trusty Tahr installation 

```shell
# dependent packages
$ sudo apt-get install cmake libgsl0-dev libboost-regex-dev libhdf5-dev libatlas-base-dev python-dev python-pip
$ sudo pip install cython

$ wget https://github.com/ecell/ecell4/archive/master.zip   
$ unzip master.zip
$ cd ecell4-master
# By default install.sh tries to install E-Cell4 into /usr/local, in this case you need to use sudo.
# In the following command, we install E-Cell4 into $HOME/ecell4. In this case you do NOT need to use sudo.
$ PREFIX=$HOME/ecell4 PYTHONPATH=/path/to/lib/python2.7/site-packages ./install.sh py2
```

## Running E-Cell4

```
# If you set PREFIX to $HOME/ecell4, make sure to append $HOME/ecell4/lib to LD_LIBRARY_PATH 
$ LD_LIBRARY_PATH=$HOME/ecell4/lib:$LD_LIBRARY_PATH PYTHONPATH=$HOME/ecell4/lib/python2.7/site-packages python
Python 2.7.6 (default, Mar 22 2014, 22:59:56) 
[GCC 4.8.2] on linux2
Type "help", "copyright", "credits" or "license" for more information.
>>> from ecell4.core import *
>>> sp = Species("B.A.C")
>>> print sp.serial()
A.B.C
>>> 
```

## Example

Here is an extremely simple example with a reversible binding reaction:


    %matplotlib inline
    import numpy
    from ecell4 import *
    
    with reaction_rules():
        A + B == C | (0.01, 0.3)
    
    y = run_simulation(
        numpy.linspace(0, 10, 100), {'A': 60, 'B': 60}, solver='ode')


![png](https://raw.githubusercontent.com/ecell/ecell4/develop/docs/output_7_0.png)


## Running E-Cell4 with IPython notebook (by docker)

We support docker images too.
If you use docker, you can easily test E-Cell4.

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

Typically, it is 192.168.59.103, then please open 192.168.59.103:8888 with your favorite browser.
(But it could get changed by Virtualbox's DHCP implementation.)

### Docker (Linux)

```shell
$ sudo docker pull ecell4/ecell4:develop
$ sudo docker run --rm -i -t -p 8888:8888 ecell/ecell4:develop
```

Open localhost:8888 with your favorite browser.


## IPython notebooks (tutorials) for E-Cell4

Please see http://nbviewer.ipython.org/github/ecell/ecell4/blob/develop/ipynb/index.ipynb
