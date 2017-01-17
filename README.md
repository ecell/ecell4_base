E-Cell System version 4
=======================

[![Binder](http://mybinder.org/badge.svg)](http://mybinder.org/repo/ecell/ecell4-notebooks)
[![Build Status](https://travis-ci.org/ecell/ecell4.svg?branch=develop)](https://travis-ci.org/ecell/ecell4)
[![Build status](https://ci.appveyor.com/api/projects/status/github/ecell/ecell4?svg=true)](https://ci.appveyor.com/project/kaizu/ecell4)
[![Documentation Status](https://readthedocs.org/projects/ecell4/badge/?version=latest)](http://ecell4.readthedocs.org/en/latest/?badge=latest)
[![PyPI](https://img.shields.io/pypi/v/ecell.svg)](https://pypi.python.org/pypi/ecell)
[![License: GPL v2](https://img.shields.io/badge/license-GPL%20v2-blue.svg)](https://github.com/ecell/ecell4/blob/master/licenses/LICENSE)
[![Slack Status](https://img.shields.io/badge/chat-on%20slack-50baa6.svg)](https://ecell-project.herokuapp.com/)
<!---[![Slack Status](https://ecell-project.herokuapp.com/badge.svg)](https://ecell-project.herokuapp.com/)--->

What is E-Cell System?
----------------------

E-Cell System is, a software platform for modeling, simulation and analysis of complex, heterogeneous and multi-scale systems like the cell.
E-Cell has multi-algorithm, multi-timescale and multi-spatial-representation as its central feature.

Quick start
===========

You can try E-Cell4 without installing it on your computer’s hard drive.

You can do this by just clicking [![Binder](http://mybinder.org/badge.svg)](http://mybinder.org/repo/ecell/ecell4-notebooks).

All you need to do is just running the cells in each of the example Jupyter notebooks.

Installation and usage
======================

- [Installation](#installation)
  - [Requirements](#requirements)
  - [Windows](#windows)
  - [Mac or Linux](#mac-or-linux)

- [How to try E-Cell4 examples](#how-to-try-e-cell4-examples)

- [Simple examples](#simple-examples)
  - [A reversible binding reaction](#a-reversible-binding-reaction)
  - [Particle tracking on a spherical surface](#particle-tracking-on-a-spherical-surface)

- [Docker image for E-Cell4](#docker-image-for-e-cell4)

Installation
------------

### Requirements

#### Minimum requirements

- Python (2.7 and 3.4, 3.5, 3.6 both major versions are supported [3.4 is only supported on Linux, 3.6 is only supported on Mac])
- pip (8.1 or later)
- hdf5 (required only on **Windows**.)

#### Optional requirements

We strongly recommend that you run E-Cell4 with [Jupyter Notebook](http://jupyter.org/).
Some E-Cell4 functions (for visualization, datastore) optionally depend on

- matplotlib (**1.5.1** or later)
- ffmpeg
- pandas

### Windows

We recommend that you install [Miniconda](http://conda.pydata.org/miniconda.html) to manage Python packages.
**Note that we do not support Python2.7 64bit for Windows.**

- Install Miniconda for Windows from http://conda.pydata.org/miniconda.html
- Run the following commands on command prompt

    ```shell
    conda install hdf5 matplotlib notebook
    pip install ecell
    ```

- (**Important**) E-Cell4 for Windows needs the latest `hdf5`. If there's any problem, please update the version of hdf5.

    ```shell
    conda update hdf5
    ```

If you use animated visualization with E-Cell4, please install [ffmpeg windows build](http://ffmpeg.zeranoe.com/builds/) and add its path to your **USER** PATH enviromental variable.

### Mac or Linux

We recommend that you install [Miniconda](http://conda.pydata.org/miniconda.html) to manage Python packages.
After installing Miniconda, run the following commands in your terminal.

(NOTICE for Mac users) We do not provide **Python3.5 whl for Mac**. Instead we provide **Python3.6 whl for Mac**. To use Python3.6 enviroment, please refer to http://conda.pydata.org/docs/py2or3.html . Continuum.io already offers Python3.6 conda packages.

```shell
# After installing Miniconda2 or Miniconda3 (Here we assume that you installed Miniconda3).
~/miniconda2/bin/conda install matplotlib notebook

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
Here we assume that you installed Miniconda(Python3.5) to C:¥Miniconda3

```shell
cd C:¥Miniconda3¥Scripts
wget https://github.com/ecell/ecell4-notebooks/archive/master.zip -OutFile master.zip
Expand-Archive master.zip
.¥jupyter-notebook.exe .¥master¥ecell4-notebooks-master¥
```

### Mac or Linux
Here we assume that you installed Miniconda(Python3.5) to ~/miniconda3

```shell
wget https://github.com/ecell/ecell4-notebooks/archive/master.zip
unzip master.zip
cd ecell4-notebooks-master
~/miniconda3/bin/jupyter-notebook
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
>>> print(sp.serial())
B.A.C
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

![png](./readthedocs/images/output_7_0.png)

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

![png](./readthedocs/images/hairball.png)

Docker image for E-Cell4
----------------------------

If you're familiar with Docker, the following commands should work in most cases:

```shell
docker pull ecell/ecell4
docker run -d -p 8888:8888 ecell/ecell4 start-notebook.sh --NotebookApp.token=''
```

and open a web browser to `http://localhost:8888` .

Our Docker image is based on **Minimal Jupyter Notebook Stack**. See https://github.com/jupyter/docker-stacks/tree/master/base-notebook or [Our Wiki page](https://github.com/ecell/ecell4/wiki/Security-in-the-Docker-Jupyter-notebook-server) for more details on the Docker command options.

Licensing terms
===============

This product is licensed under the terms of the [GNU General Public License v2](https://github.com/ecell/ecell4/blob/master/licenses/LICENSE),
See [NOTICE](https://github.com/ecell/ecell4/blob/master/licenses/NOTICE.txt) for the software included in this product.

- Copyright (c) 2010-, RIKEN

All rights reserved.
