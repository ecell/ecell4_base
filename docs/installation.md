Installation
=================

Quick start with docker
- [Docker container for E-Cell System version4](#docker-container-for-e-cell-system-version4)
  - [Windows or Mac](#windows-or-mac)
  - [Linux](#linux)

Installation and usage
- [Requirements](#requirements)
  - [Minimum requirements](#minimum-requirements)
  - [Optional requirements](#optional-requirements)
  - [Build requirements](#build-requirements)

- [Windows](#windows)
  - [Python2 series](#python2-series)
  - [Python3 series](#python3-series)
- [Mac](#mac)
  - [Homebrew users](#homebrew-users)
  - [Pip users](#pip-users)
- [Linux](#Linux)
  - [Custom shell script for Ubuntu](#custom-shell-script-for-ubuntu)
  - [Linuxbrew](#linuxbrew)

- [Simple examples](#simple-examples)


Docker container for E-Cell System version4
-------------------------------------------

If you have docker environment, you can easily try E-Cell4.
You can pull E-Cell4 container with `docker pull ecell/ecell4`.

After the following steps, you should see [Jupyter Notebook](http://jupyter.org/) up and running (and E-Cell4 tutorials) in your web browser.

### Windows or Mac

1. Install [Docker Toolbox](https://www.docker.com/toolbox).
2. Run **Docker Quickstart Terminal**.
3. Run `docker run -d -p 443:8888 ecell/ecell4` in your command prompt or terminal.
4. Open **192.168.99.100:443** with your favorite web browser.

### Linux

1. Install docker.
2. Run the following commands in your terminal.

    ```shell
    sudo docker pull ecell/ecell4
    sudo docker run -d -p 443:8888 ecell/ecell4
    ```

3. Open **localhost:443** with your favorite web browser.


Installation and usage
======================
Requirements
------------
### Minimum requirements
- Python 2.7 series or Python 3.5 series
- pip

### Optional requirements
We strongly recommend that you run E-Cell4 from [Jupyter Notebook](http://jupyter.org/).
And some E-Cell4 functions (for datastore, visualization) optionaly depend on
- hdf5
- matplotlib **1.5.1** and later
- ffmpeg or avconv
- pandas

### Build requirements
If you build E-Cell4 from source code, you need to install these software.
- cmake
- boost
- gsl
- hdf5

Windows
-------

Please use 32bit Python, even if you use 64bit Windows.
We have NOT supported 64bit Python yet.

### Python2 series

- [Python 2.7.11(**32bit**)](https://www.python.org/ftp/python/2.7.11/python-2.7.11.msi)
- [HDF5-1.8.16(**32-bit built with VS2012**)](http://www.hdfgroup.org/ftp/HDF5/current/bin/windows/extra/hdf5-1.8.16-win32-vs2012-shared.zip)
- [Visual C++ Compiler for Python 2.7](http://aka.ms/vcpython27)

Please add python.exe, pip.exe path and `C:\Program Files (x86)\HDF_Group\HDF5\1.8.16\bin` to your **USER** PATH enviromental variable.

And run the following commands with command prompt.
matplotlib depends on numpy. It takes some time to build numpy, please be patient.
Although jupyter is optional, we strongly recommend that you run E-Cell4 with jupyter.

```shell
pip install https://ci.appveyor.com/api/buildjobs/59qrnnjpqgwdrot5/artifacts/python/dist/ecell4-4.0.0b2-cp27-none-win32.whl
pip install -U matplotlib
pip install -U jupyter
```

If you use animated visualization for E-Cell4, please install [ffmpeg windows build](http://ffmpeg.zeranoe.com/builds/) and add its path to your **USER** PATH enviromental variable too.

### Python3 series

- [Python 3.5.1(**32bit**)](https://www.python.org/ftp/python/3.5.1/python-3.5.1.msi)
- [HDF5-1.8.16(**32-bit built with VS2015**)](http://www.hdfgroup.org/ftp/HDF5/current/bin/windows/extra/hdf5-1.8.16-win32-vs2015-shared.zip)

Please add python.exe, pip.exe path and `C:\Program Files (x86)\HDF_Group\HDF5\1.8.16\bin` to your **USER** PATH enviromental variable.
Next download numpy-1.10.4+vanilla-cp35-none-win32.whl and matplotlib-1.5.0-cp35-none-win32.whl from http://www.lfd.uci.edu/~gohlke/pythonlibs/
and run the following commands with command prompt.

```shell
pip install https://ci.appveyor.com/api/buildjobs/jpyueyasgwsannch/artifacts/python/dist/ecell4-4.0.0b2-cp35-none-win32.whl
pip install numpy-1.10.4+vanilla-cp35-none-win32.whl
pip install matplotlib-1.5.1-cp35-none-win32.whl
pip install -U jupyter
```

If you use animated visualization for E-Cell4, please install [ffmpeg windows build](http://ffmpeg.zeranoe.com/builds/) and add its path to your **USER** PATH enviromental variable too.

Mac
---

### homebrew users
First, we recommend that Mac users start E-Cell4 with [homebrew-ecell4](https://github.com/ecell/homebrew-ecell4).
This homebrew formula includes ffmpeg (for animated visualization).
Please run the following commands in your terminal.

```shell
brew tap ecell/ecell4
brew install ecell4
# Mac default matplotlib is too old for E-Cell4, you need to update it with the following options.
pip install -U matplotlib --user
pip install -U jupyter --user
# path config for homebrew-ecell4
mkdir -p ~/Library/Python/2.7/lib/python/site-packages
echo '/usr/local/lib/python2.7/site-packages' >> ~/Library/Python/2.7/lib/python/site-packages/homebrew.pth
# path config for --user installed Python packages
echo 'export PYTHONPATH=~/Library/Python/2.7/lib/python/site-packages:$PYTHONPATH' >> ~/.bashrc
echo 'export PATH=~/Library/Python/2.7/bin:$PATH' >> ~/.bashrc
source ~/.bashrc
```

### pip users
We also have Python wheel files for E-Cell4.
But the wheel distribution does NOT include ffmpeg.

```shell
# please select appropriate whl file for your Python version
pip install THEWHEELURL.whl --user
# Mac default matplotlib is too old for E-Cell4, you need to update it with the following options.
pip install -U matplotlib --user
pip install -U jupyter --user
# path config for --user installed Python packages
echo 'export PYTHONPATH=~/Library/Python/2.7/lib/python/site-packages:$PYTHONPATH' >> ~/.bashrc
echo 'export PATH=~/Library/Python/2.7/bin:$PATH' >> ~/.bashrc
source ~/.bashrc
```

Linux
-----

### Custom shell script for Ubuntu
First, we recommend that Ubuntu(14.04 LTS) users start E-Cell4 with our custom shell script named **install.sh**.
Please run the following commands in your terminal.
(In these commands, we assume that you are **root* user.)

```shell
apt-get install -y software-properties-common libav-tools python-dev libfreetype6-dev libpng-dev pkg-config pandoc wget cmake g++ libboost-dev libgsl0-dev libhdf5-serial-dev libboost-regex-dev python python-numpy python-scipy python-pip python-zmq
add-apt-repository ppa:mc3man/trusty-media -y; apt-get update; apt-get install -y ffmpeg

pip install cython jupyter matplotlib

wget THE_RELEASE_URL.tar.gz
tar xf THE_RELEASE_URL.tar.gz
cd ecell4; export PREFIX=/usr/local; export PYTHONPATH=/usr/local/lib/python2.7/site-packages:$PYTHONPATH; ./install.sh --python2 --hdf5
```

### Linuxbrew
[E-Cell4 homebrew formula](https://github.com/ecell/homebrew-ecell4) also can be used for [Linuxbrew](http://linuxbrew.sh/).
If you do NOT use Ubuntu, please try Linuxbrew instead of **install.sh**.

Simple examples
---------------

Here are two extremely simple examples, See http://ecell.github.io/ecell4/ for more details on running E-Cell4.

```
Python 2.7.6 (default, Mar 22 2014, 22:59:56)
[GCC 4.8.2] on linux2
Type "help", "copyright", "credits" or "license" for more information.
>>> from ecell4.core import *
>>> sp = Species("B.A.C")
>>> print sp.serial()
A.B.C
>>>
```

### A reversible binding reaction

```python
%matplotlib inline
import numpy
from ecell4 import *

with reaction_rules():
    A + B == C | (0.01, 0.3)

y = run_simulation(
    numpy.linspace(0, 10, 100), {'A': 60, 'B': 60}, solver='ode')
```

![png](https://raw.githubusercontent.com/ecell/ecell4/master/docs/images/output_7_0.png)

### Particle tracking on a spherical surface

```python
%matplotlib inline
from ecell4 import *

with species_attributes():
    A | {'D': '1', 'location': 'M'}

surface = Sphere(Real3(0.5, 0.5, 0.5), 0.5).surface()
obs = FixedIntervalTrajectoryObserver(1e-4)
run_simulation(
    0.4, y0={'A': 10}, structures={'M': surface},
    solver='spatiocyte', observers=obs, return_type=None)

viz.plot_trajectory(obs, interactive=False)
```

![png](https://raw.githubusercontent.com/ecell/ecell4/master/docs/images/hairball.png)
