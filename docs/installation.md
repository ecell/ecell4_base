- [Installation](#installation)
    - [Windows](#windows-installation)
    - [Mac OS X](#mac-os-x-installation)
    - [Ubuntu](#ubuntu-linux-installation)
    - [CentOS](#centos-linux-installation)
- [Running E-Cell4](#running-e-cell4)
- [Dockerized E-Cell4 Jupyter notebooks](#dockerized-e-cell4-jupyter-notebooks)
    - [For Windows and Mac](#for-windows-and-mac)
    - [For Linux](#for-linux)

Installation
------------

### Windows installation

#### Requirements

Please use 32bit Python, even if you use 64bit Windows.
We don't support 64bit Python

- Python 2.7.10(**32bit**) https://www.python.org/ftp/python/2.7.10/python-2.7.10.msi
- HDF5-1.8.16 Pre-built Binary(**32-bit**) https://www.hdfgroup.org/ftp/HDF5/current/bin/windows/extra/hdf5-1.8.16-win32-vs2012-shared.zip

Please add `C:\Python27`, `C:\Python27\Scripts` and `C:\Program Files (x86)\HDF_Group\HDF5\1.8.16\bin` to your **PATH** enviromental variable.

And run following command with command prompt.
```
pip install https://github.com/ecell/ecell4/releases/download/4.0.0-beta2/ecell4-4.0.0b2-cp27-none-win32.whl
```

#### Jupyter for Windows
We recommend you run E-Cell4 models from Jupyter notebook.
Below is Jupyter notebook(and matplotlib) installation for Windows.

- Install [Visual C++ Compiler for Python 2.7](http://www.microsoft.com/en-us/download/details.aspx?id=44266)
- Install Jupyter notebook and matplotlib

  ```
  pip install -U jupyter
  pip install matplotlib
  ```

matplotlib depends on numpy. It takes some time to build numpy, please be patient.

### Mac OS X installation

```shell
ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
brew tap ecell/ecell4
brew install ecell4 --HEAD
```

#### Jupyter for Mac
We recommend you run E-Cell4 models from Jupyter notebook.
Below is Jupyter notebook(and matplotlib) installation for Mac.

```shell
sudo python get-pip.py
sudo pip install -U matplotlib
sudo pip install -U jupyter
```

### Ubuntu Linux installation

```shell
$ sudo apt-get install build-essential curl git m4 ruby texinfo libbz2-dev libcurl4-openssl-dev libexpat-dev libncurses-dev zlib1g-dev python-dev python-setuptools
$ ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/linuxbrew/go/install)"
$ echo 'export PATH="$HOME/.linuxbrew/bin:$PATH"' >> ~/.bash_profile
$ source ~/.bash_profile
$ brew tap ecell/ecell4
$ brew install ecell4 --HEAD
```

### CentOS Linux installation

```shell
$ sudo yum groupinstall 'Development Tools'
$ sudo yum install curl git m4 ruby texinfo bzip2-devel curl-devel expat-devel ncurses-devel zlib-devel python-devel python-setuptools
$ ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/linuxbrew/go/install)"
$ echo 'export PATH="$HOME/.linuxbrew/bin:$PATH"' >> ~/.bash_profile
$ source ~/.bash_profile
$ ln -s `which gcc` `brew --prefix`/bin/gcc-4.8
$ ln -s `which g++` `brew --prefix`/bin/g++-4.8
$ brew tap ecell/ecell4
$ brew install ecell4 --HEAD
```

Running E-Cell4
---------------

### How to run Python for E-Cell4 (Windows)
1. run command prompt
2. run Python
  ```
  python
  ```

### How to run Python for E-Cell4 (Mac)
1. run terminal
2. run Python with PYTHONPATH
  ```
  PYTHONPATH=/usr/local/Cellar/ecell4/HEAD/libexec/lib/python2.7/site-packages python
  ```

### How to run Python for E-Cell4 (Ubuntu)

If you set PREFIX to $HOME/ecell4, run Python with
```
LD_LIBRARY_PATH=$HOME/ecell4/lib:$LD_LIBRARY_PATH PYTHONPATH=$HOME/ecell4/lib/python2.7/site-packages python
```
in case with Python3
```
LD_LIBRARY_PATH=$HOME/ecell4/lib:$LD_LIBRARY_PATH PYTHONPATH=$HOME/ecell4/lib/python3.4/site-packages python3
```

### Simple examples

Here are two extremely simple examples, See http://ecell4.readthedocs.org/en/latest/tutorials/ for more details on running E-Cell4.

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

![png](https://raw.githubusercontent.com/ecell/ecell4/master/docs/output_7_0.png)


Dockerized E-Cell4 Jupyter notebooks
------------------------------------

If you use docker, you can easily try E-Cell4.
You can pull E-Cell4 container with `docker pull ecell/ecell4`

### For Windows and Mac

1. Install [Docker Toolbox](https://www.docker.com/toolbox)
2. Run Kitematic
3. Search with **ecell4**, and create ecell4 container

  ![png](https://raw.githubusercontent.com/ecell/ecell4/master/docs/kitematic1.png)

4. Open the **ACCESS URL** in **IP & PORTS** with your web browser

  ![png](https://raw.githubusercontent.com/ecell/ecell4/master/docs/kitematic2.png)

### For Linux

```shell
$ sudo docker pull ecell/ecell4
$ sudo docker run -d -p 443:8888 ecell/ecell4
```

You'll now be able to E-Cell4 notebooks at https://THE_IP_RUNNING_DOCKER:443
