E-Cell System version 4
=======================

[![Build Status](https://travis-ci.org/ecell/ecell4.svg?branch=develop)](https://travis-ci.org/ecell/ecell4)
[![Build status](https://ci.appveyor.com/api/projects/status/github/ecell/ecell4?svg=true)](https://ci.appveyor.com/project/kaizu/ecell4)
[![Documentation Status](https://readthedocs.org/projects/ecell4/badge/?version=latest)](http://ecell4.readthedocs.org/en/latest/?badge=latest)
[![GitHub release](https://img.shields.io/github/release/ecell/ecell4.svg)](https://github.com/ecell/ecell4/releases)
[![GitHub license](https://img.shields.io/github/license/ecell/ecell4.svg)](https://github.com/ecell/ecell4/LICNESE)
[![Gitter](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/ecell/ecell4?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge)

What is E-Cell System?
----------------------

E-Cell System is, a software platform for modeling, simulation and analysis of complex, heterogeneous and multi-scale systems like the cell.

Installation and usage
======================

- [Docker container for E-Cell System version4](#docker-container-for-e-cell-system-version4)
  - [Windows or Mac](#windows-or-mac)
  - [Linux](#linux)

- [Installation](#installation)
  - [Requirements](#requirements)
  - [Windows](#windows)
  - [Mac](#mac)
  - [Linux](#linux-1)

- [How to open E-Cell4 Jupyter notebooks](#how-to-open-e-cell4-jupyter-notebooks)

- [Simple examples](#simple-examples)
  - [A reversible binding reaction](#a-reversible-binding-reaction)
  - [Particle tracking on a spherical surface](#particle-tracking-on-a-spherical-surface)

Docker container for E-Cell System version4
-------------------------------------------

If you have docker environment, you can easily try E-Cell4.
You can pull E-Cell4 container with `docker pull ecell/ecell4`.

After the following steps, you should see [Jupyter Notebook](http://jupyter.org/) up and running (and E-Cell4 tutorials) in your web browser.

### Windows or Mac

1. Install [Docker Toolbox](https://www.docker.com/toolbox).
2. Run **Docker Quickstart Terminal**.
3. Run the following commands

    ```shell
    docker pull ecell/ecell4
    docker run -dp 443:8888 ecell/ecell4
    ```

4. Open **192.168.99.100:443** with your web browser.

### Linux

1. Install Docker.
2. Run the following commands in your terminal.

    ```shell
    sudo docker pull ecell/ecell4
    sudo docker run -dp 443:8888 ecell/ecell4
    ```

3. Open **localhost:443** with your web browser.


Installation
------------

### Requirements

#### Minimum requirements
- Python or **32bit** Miniconda for Windows (2.7, 3.* both versions are supported)
- pip
- hdf5 (E-Cell4 **for Windows** works only for version 1.8.16)

#### Optional requirements
We strongly recommend that you run E-Cell4 with [Jupyter Notebook](http://jupyter.org/).
And some E-Cell4 functions (for visualization, datastore) optionaly depend on
- matplotlib (**1.5.1** and later)
- ffmpeg
- pandas

### Windows

Please use **32bit** [Miniconda](http://conda.pydata.org/miniconda.html), even if you use 64bit Windows.
We have NOT supported 64bit Python yet.
Python 2.7, 3.5 both are supported.

- Install **32bit** Miniconda for Windows from http://conda.pydata.org/miniconda.html
- Run the follwing commands on command prompt (if you use Python3.5, please replace the target of ```pip install``` to the whl for 3.5)
- (**Important**) E-Cell4 for Windows works only for hdf5 version **1.8.16**. Please check the version of hdf5, even if you installed hdf5 before with conda.

    ```shell
    conda install hdf5 notebook matplotlib
    pip install https://github.com/ecell/ecell4/releases/download/4.0.0/ecell-4.0.0-cp27-none-win32.whl
    ```

Although Jupyter Notebook is optional, we strongly recommend that you run E-Cell4 with jupyter.
If you use animated visualization for E-Cell4, please install [ffmpeg windows build](http://ffmpeg.zeranoe.com/builds/) and add its path to your **USER** PATH enviromental variable.

### Mac

Please run the following commands in your terminal.

```shell
# Please download E-Cell4 whl file for your Python version from https://github.com/ecell/ecell4/releases , here we downloaded a whl for Python27
pip install  --user ecell-4.0.0-cp27-none-macosx_10_6_intel.macosx_10_9_intel.macosx_10_9_x86_64.macosx_10_10_intel.macosx_10_10_x86_64.whl

# Mac default matplotlib is too old for E-Cell4, you need to update it with the following options.
pip install -U --user matplotlib
pip install -U --user jupyter

# path config for --user installed Python packages
echo 'export PYTHONPATH=~/Library/Python/2.7/lib/python/site-packages:$PYTHONPATH' >> ~/.bashrc
echo 'export PATH=~/Library/Python/2.7/bin:$PATH' >> ~/.bashrc
source ~/.bashrc

# If you use animation support. (Install ffmpeg with homebrew)
brew install ffmpeg
```


### Linux

Please run the following commands with root privilege.

```shell
apt-get install libgsl0-dev libhdf5-dev wget
wget https://bootstrap.pypa.io/get-pip.py
# If you use Python3, replace python to python3
python get-pip.py
# If you use Python3.*, replace the whl for Python3.*
pip install https://github.com/ecell/ecell4/releases/download/4.0.0/ecell-4.0.0-cp27-cp27mu-manylinux1_x86_64.whl

# The latest matplotlib and jupyter. If you use Python3, replace those for Python3.
apt-get install python-dev libfreetype6-dev libpng-dev pkg-config python-numpy pandoc
pip install matplotlib jupyter

# Optional requirement (animation visualization)
apt-get install libav-tools
```


How to open E-Cell4 Jupyter notebooks
-------------------------------------

### Windows
Please replace the CONDA_INSTALL_FOLDER with the folder you installed Miniconda.
For example **C:¥Miniconda27**.

```shell
cd the CONDA_INSTALL_FOLDER
cd ecell4ipynb
jupyter-notebook
```

### Mac

```shell
### in the case of Python27
~/Library/Python/2.7/ecell4ipynb
jupyter-notebook
```

### Linux
```shell
cd /usr/local/ecell4ipynb
jupyter-notebook
```

Simple examples
---------------

Here are two extremely simple examples, See http://ecell4.readthedocs.org for more details on running E-Cell4.

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

surface = Sphere(ones() * 0.5, 0.5).surface()
obs = FixedIntervalTrajectoryObserver(1e-4)
run_simulation(
    0.4, y0={'A': 10}, structures={'M': surface},
    solver='spatiocyte', observers=obs, return_type=None)

viz.plot_trajectory(obs, interactive=False)
```

![png](https://raw.githubusercontent.com/ecell/ecell4/master/docs/images/hairball.png)

Licensing terms
===============

This project is licensed under the terms of the GNU General Public License v2.
See [LICENSE](https://github.com/ecell/ecell4/blob/master/LICENSE) for the project license.

- Copyright (c) 2010-, RIKEN

All rights reserved.
