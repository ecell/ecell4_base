E-Cell System version 4
=======================

[![Binder](http://mybinder.org/badge.svg)](http://mybinder.org/repo/ecell/ecell4-notebooks)
[![Build Status](https://travis-ci.org/ecell/ecell4.svg?branch=develop)](https://travis-ci.org/ecell/ecell4)
[![Build status](https://ci.appveyor.com/api/projects/status/github/ecell/ecell4?svg=true)](https://ci.appveyor.com/project/kaizu/ecell4)
[![Documentation Status](https://readthedocs.org/projects/ecell4/badge/?version=latest)](http://ecell4.readthedocs.org/en/latest/?badge=latest)
[![PyPI](https://img.shields.io/pypi/v/ecell.svg)](https://pypi.python.org/pypi/ecell)
[![License: GPL v2](https://img.shields.io/badge/license-GPL%20v2-blue.svg)](https://github.com/ecell/ecell4/blob/master/licenses/LICENSE)
[![Slack Status](https://ecell-project.herokuapp.com/badge.svg)](https://ecell-project.herokuapp.com/)

What is E-Cell System?
----------------------

E-Cell System is, a software platform for modeling, simulation and analysis of complex, heterogeneous and multi-scale systems like the cell.

Installation and usage
======================

- [Docker container for E-Cell4](#docker-container-for-e-cell4)

- [Installation](#installation)
  - [Requirements](#requirements)
  - [Windows](#windows)
  - [Mac Linux](#mac-linux)

- [How to try E-Cell4 examples](#how-to-try-e-cell4-examples)

- [Simple examples](#simple-examples)
  - [A reversible binding reaction](#a-reversible-binding-reaction)
  - [Particle tracking on a spherical surface](#particle-tracking-on-a-spherical-surface)

Docker container for E-Cell4
----------------------------

If you have docker environment, you can easily try E-Cell4.
You can pull E-Cell4 container with `docker pull ecell/ecell4`.

After the following steps, you should see [Jupyter Notebook](http://jupyter.org/) up and running (and E-Cell4 tutorials) in your web browser.

1. Install [Docker](https://www.docker.com/products/docker).
2. Run Docker.
3. Run the following commands

    ```shell
    docker pull ecell/ecell4
    docker run -dP ecell/ecell4
    ```

4. Check which port is used by E-Cell4 docker with `docker ps` command.

    ```shell
    docker ps
    CONTAINER ID        IMAGE               COMMAND                  CREATED             STATUS              PORTS                     NAMES
    82b90db240f5        ecell/ecell4        "/bin/sh -c 'jupyter-"   7 seconds ago       Up 6 seconds        0.0.0.0:32768->8888/tcp   clever_raman
    ```

5. Open the PORT in `docker ps` command with your web browser. In the case of the above example, you will open `0.0.0.0:32768`


Installation
------------

### Requirements

#### Minimum requirements
- Python or **32bit** Miniconda for Windows (2.7 and 3.4, 3.5 both major versions are supported [3.4 is only supported on Linux])
- pip (8.1 and later)
- hdf5 (required only on **Windows**. works only for **version 1.8.16**)

#### Optional requirements
We strongly recommend that you run E-Cell4 with [Jupyter Notebook](http://jupyter.org/).
And some E-Cell4 functions (for visualization, datastore) optionaly depend on
- matplotlib (**1.5.1** and later)
- ffmpeg
- pandas

### Windows

Please use **32bit** [Miniconda](http://conda.pydata.org/miniconda.html), even if you use 64bit Windows.
We have NOT supported 64bit Python yet.

- Install **32bit** Miniconda for Windows from http://conda.pydata.org/miniconda.html
- Run the follwing commands on command prompt
- (**Important**) E-Cell4 for Windows works only for hdf5 version **1.8.16**. Please check the version of hdf5, even if you installed hdf5 before with conda.

    ```shell
    conda install hdf5 notebook matplotlib
    pip install ecell
    ```

Although Jupyter Notebook is optional, we strongly recommend that you run E-Cell4 with jupyter.
If you use animated visualization with E-Cell4, please install [ffmpeg windows build](http://ffmpeg.zeranoe.com/builds/) and add its path to your **USER** PATH enviromental variable.

### Mac Linux

We recommend that you install [Miniconda](http://conda.pydata.org/miniconda.html) to manage Python packages.
After installing Miniconda, run the following commands in your terminal.

```shell
# After installing Miniconda2 or Miniconda3 (Here we assume that you installed Miniconda2).
~/miniconda2/bin/conda install matplotlib jupyter

# Download E-Cell4 whl file for your Python version from https://github.com/ecell/ecell4/releases before running this command.
~/miniconda2/bin/pip install ecell

# If you want animation support, install ffmpeg with homebrew
brew install ffmpeg
# or if you use Ubuntu Linux
# apt install ffmpeg
```

How to try E-Cell4 examples
---------------------------
Here we download example notebooks from https://github.com/ecell/ecell4-notebooks and open it with Jupyter Notebook.

### Windows
Open powershell and run these commands.
Here we assume that you installed Miniconda(Python2.7) to C:¥Miniconda2

```shell
cd C:¥Miniconda2¥Scripts
wget https://github.com/ecell/ecell4-notebooks/archive/master.zip -OutFile master.zip
Expand-Archive master.zip
.¥jupyter-notebook.exe .¥master¥ecell4-notebooks-master¥
```

### Mac Linux
Here we assume that you installed Miniconda(Python2.7) to ~/miniconda2

```shell
wget https://github.com/ecell/ecell4-notebooks/archive/master.zip
unzip master.zip
cd ecell4-notebooks-master
~/miniconda2/bin/jupyter-notebook
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

This product is licensed under the terms of the [GNU General Public License v2](https://github.com/ecell/ecell4/blob/master/licenses/LICENSE),
See [NOTICE](https://github.com/ecell/ecell4/blob/master/licenses/NOTICE.txt) for the software included in this product.

- Copyright (c) 2010-, RIKEN

All rights reserved.
