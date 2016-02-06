Installation and usage
======================

- [Docker container for E-Cell System version4](#docker-container-for-e-cell-system-version4)
  - [Windows or Mac](#windows-or-mac)
  - [Linux](#linux)

- [Installation](#installation)
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
    - [Ubuntu Linux](#ubuntu-linux)
    - [Linuxbrew](#linuxbrew)

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
3. Run `docker run -d -p 443:8888 ecell/ecell4` in the terminal.
4. Open **192.168.99.100:443** with your web browser.

### Linux

1. Install Docker.
2. Run the following commands in your terminal.

    ```shell
    sudo docker pull ecell/ecell4
    sudo docker run -d -p 443:8888 ecell/ecell4
    ```

3. Open **localhost:443** with your web browser.


Installation
------------

### Requirements

#### Minimum requirements
- Python 2.7 or 3.4(on Linux) 3.5(on Windows, Mac) 
- pip

#### Optional requirements
We strongly recommend that you run E-Cell4 with [Jupyter Notebook](http://jupyter.org/).
And some E-Cell4 functions (for visualization, datastore) optionaly depend on
- matplotlib (**1.5.1** and later)
- ffmpeg or avconv
- hdf5
- pandas

### Build requirements
If you build E-Cell4 from source code, you need to install these software.
- boost (**1.59** and earlier)
- cmake
- gsl
- hdf5

### Windows

Please use 32bit Python2.7 or 3.5, even if you use 64bit Windows.
We have NOT supported 64bit Python yet.

- Install **32bit** Miniconda for Windows from http://conda.pydata.org/miniconda.html
- Run the follwing commands on command prompt (if you use Python3.5, please replace the target of ```pip install``` to the whl for 3.5)

    ```shell
    conda update pip
    conda install hdf5 jupyter matplotlib
    pip install https://github.com/ecell/ecell4/releases/download/4.0.0/ecell-4.0.0-cp27-none-win32.whl
    ``` 

Although jupyter is optional, we strongly recommend that you run E-Cell4 with jupyter.
If you use animated visualization for E-Cell4, please install [ffmpeg windows build](http://ffmpeg.zeranoe.com/builds/) and add its path to your **USER** PATH enviromental variable.

### Mac

#### homebrew users
We have homebrew formula for E-Cell4 [homebrew-ecell4](https://github.com/ecell/homebrew-ecell4).
This homebrew formula includes ffmpeg (for animated visualization).
Please run the following commands in your terminal.

```shell
brew tap ecell/ecell4
brew install ecell4

# Mac default matplotlib is too old for E-Cell4, you need to update it with the following options.
pip install -U --user matplotlib
pip install -U --user jupyter

# path config for homebrew-ecell4
mkdir -p ~/Library/Python/2.7/lib/python/site-packages
echo 'import site; site.addsitedir("/usr/local/lib/python2.7/site-packages")' >> ~/Library/Python/2.7/lib/python/site-packages/homebrew.pth

# path config for --user installed Python packages
echo 'export PYTHONPATH=~/Library/Python/2.7/lib/python/site-packages:$PYTHONPATH' >> ~/.bashrc
echo 'export PATH=~/Library/Python/2.7/bin:$PATH' >> ~/.bashrc
source ~/.bashrc
```

#### pip users
We also have Python wheel files for E-Cell4.
But the wheel distribution does NOT include ffmpeg.
If you need animation visualization, please install it yourself. 

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
```

### Linux

#### Ubuntu Linux
We have tested the release files on Ubuntu Linux 14.04 and 15.10.
If you use Ubuntu please run the following commands.

```shell
# If you use Python3 please replace python-pip to python3-pip
sudo apt-get install libboost-dev libgsl0-dev libhdf5-dev python-pip
# If you use Python3 please replace the whl for Python3
pip install --user https://github.com/ecell/ecell4/releases/download/4.0.0/ecell-4.0.0-cp27-none-linux_x86_64.whl

# The latest matplotlib and jupyter
sudo apt-get install libfreetype6-dev libpng-dev pkg-config python-numpy pandoc
pip install --user matplotlib jupyter

# Optional requirement (animation visualization) for 14.04
sudo apt-get install libav-tools
# If you use 15.10 you can use ffmpeg deb package
#sudo apt-get install -y ffmpeg

# path config for --user installed Python packages
echo 'export PYTHONPATH=~/.local/lib/python2.7/site-packages:$PYTHONPATH' >> ~/.bashrc
echo 'export PATH=~/.local/bin:$PATH' >> ~/.bashrc
source ~/.bashrc
```

#### Linuxbrew

[E-Cell4 homebrew formula](https://github.com/ecell/homebrew-ecell4) also can be used for [Linuxbrew](http://linuxbrew.sh/).
If you do NOT use Ubuntu, please try Linuxbrew instead.

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

surface = Sphere(Real3(0.5, 0.5, 0.5), 0.5).surface()
obs = FixedIntervalTrajectoryObserver(1e-4)
run_simulation(
    0.4, y0={'A': 10}, structures={'M': surface},
    solver='spatiocyte', observers=obs, return_type=None)

viz.plot_trajectory(obs, interactive=False)
```

![png](https://raw.githubusercontent.com/ecell/ecell4/master/docs/images/hairball.png)
