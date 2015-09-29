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
- [Running E-Cell4](#running-e-cell4)
- [Dockerized E-Cell4 IPython notebooks](#dockerized-e-cell4-ipython-notebooks)
    - [For Windows and Mac](#for-windows-and-mac)
    - [For Linux](#for-linux)

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

```shell
ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
brew tap ecell/ecell4
brew install ecell4
```

#### IPython notebook
We recommend you run E-Cell4 models from IPython notebook.
Below is IPython notebook(and matplotlib) installation for Mac.

```shell
sudo python get-pip.py
sudo pip install -U matplotlib
sudo pip install -U "ipython[notebook]"
```

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
$ sudo apt-get install cmake libgsl0-dev libboost-regex-dev libhdf5-dev cython3

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

Running E-Cell4
---------------

### Simple examples

Here are two extremely simple examples, See http://ecell4.readthedocs.org/en/develop/tutorials/ for more details on running E-Cell4.

```
# If you set PREFIX to $HOME/ecell4, make sure to append $HOME/ecell4/lib to LD_LIBRARY_PATH 
$ LD_LIBRARY_PATH=$HOME/ecell4/lib:$LD_LIBRARY_PATH PYTHONPATH=$HOME/ecell4/lib/python2.7/site-packages python
# in case with Python3
# LD_LIBRARY_PATH=$HOME/ecell4/lib:$LD_LIBRARY_PATH PYTHONPATH=$HOME/ecell4/lib/python3.4/site-packages python3

Python 2.7.6 (default, Mar 22 2014, 22:59:56) 
[GCC 4.8.2] on linux2
Type "help", "copyright", "credits" or "license" for more information.
>>> from ecell4.core import *
>>> sp = Species("B.A.C")
>>> print sp.serial()
A.B.C
>>> 
```

#### A reversible binding reaction

```python
%matplotlib inline
import numpy
from ecell4 import *

with reaction_rules():
    A + B == C | (0.01, 0.3)

y = run_simulation(
    numpy.linspace(0, 10, 100), {'A': 60, 'B': 60}, solver='ode')
```

![png](https://raw.githubusercontent.com/ecell/ecell4/develop/docs/output_7_0.png)


Dockerized E-Cell4 IPython notebooks
------------------------------------

If you use docker, you can easily try E-Cell4.
You can pull E-Cell4 container with `docker pull ecell/ecell4`

### For Windows and Mac

1. Install [Docker Toolbox](https://www.docker.com/toolbox)
2. Run Kitematic
3. Search with **ecell4**, and create ecell4 container

  ![png](https://raw.githubusercontent.com/ecell/ecell4/develop/docs/kitematic1.png)

4. Open the **ACCESS URL** in **IP & PORTS** with your web browser 

  ![png](https://raw.githubusercontent.com/ecell/ecell4/develop/docs/kitematic2.png)

### For Linux

```shell
$ sudo docker pull ecell/ecell4
$ sudo docker run -d -p 443:8888 ecell/ecell4
```

You'll now be able to E-Cell4 notebooks at https://THE_IP_RUNNING_DOCKER:443

