E-Cell System version 4
=======================

[![Build Status](https://travis-ci.org/ecell/ecell4.svg?branch=develop)](https://travis-ci.org/ecell/ecell4)
[![Gitter](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/ecell/ecell4?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge)

## What is E-Cell System?

E-Cell System, a software platform for modeling, simulation and analysis of complex, heterogeneous and multi-scale systems like the cell.

- [Docker users](#docker-users)
    - [For Windows or Mac](#for-windows-or-mac)
    - [For Linux](#for-linux)

- [Pip users](#pip-users)
    - [Windows](#windows)
    - [Mac](#mac)

- [Homebrew or Linuxbrew users](#homebrew-or-linuxbrew-users)

- [Simple examples](#simple-examples)

Docker users
------------

If you have docker environment, you can try E-Cell4 easily.
You can pull E-Cell4 container with `docker pull ecell/ecell4`.

### For Windows or Mac

1. Install [Docker Toolbox](https://www.docker.com/toolbox).
2. Run **Docker Quickstart Terminal**.
3. Run `docker run -d -p 443:8888 ecell/ecell4` in the terminal.
4. Open **192.168.99.100:443** with your favorite web browser.  
5. You should see Jupyter Notebook up and running (and E-Cell4 tutorials) in your web browser.

### For Linux

1. Install docker.
2. Run the following command.

    ```shell
    $ sudo docker pull ecell/ecell4
    $ sudo docker run -d -p 443:8888 ecell/ecell4
    ```

3. Open **localhost:443** with your favorite web browser.


Pip users
---------

### Windows

#### Requirements and installation

Please use 32bit Python, even if you use 64bit Windows.
We have NOT supported 64bit Python yet.

- [Python 2.7.11(**32bit**)](https://www.python.org/ftp/python/2.7.11/python-2.7.11.msi)
- HDF5-1.8.16 Pre-built Binary(**32-bit**) http://www.hdfgroup.org/ftp/HDF5/current/bin/windows/extra/hdf5-1.8.16-win32-vs2015-shared.zip

Please add `C:\Python27`, `C:\Python27\Scripts` and `C:\Program Files (x86)\HDF_Group\HDF5\1.8.16\bin` to your **PATH** enviromental variable.

And run following command with command prompt.
```shell
pip install https://ci.appveyor.com/api/buildjobs/aju5rykh88bb88ns/artifacts/python/dist/ecell4-4.0.0b2-cp27-none-win32.whl
```

#### Jupyter and matplotlib for Windows
We recommend you run E-Cell4 models from Jupyter notebook.
Below is Jupyter notebook(and matplotlib) installation for Windows.

- Install [Visual C++ Compiler for Python 2.7](http://aka.ms/vcpython27)
- Install Jupyter notebook and matplotlib

  ```shell
  pip install -U jupyter
  pip install matplotlib
  ```

matplotlib depends on numpy. It takes some time to build numpy, please be patient.

### Mac

#### Installation

```shell
pip install https://bintray.com/artifact/download/kozo2/generic/dist/ecell4-4.0.0b2-cp27-none-macosx_10_6_intel.macosx_10_9_intel.macosx_10_9_x86_64.macosx_10_10_intel.macosx_10_10_x86_64.whl
```

#### Jupyter and matplotlib for Mac

```shell
sudo python get-pip.py
# you can NOT install latest matplotlib into System directory, so you need to install it into USER directory
pip install -U matplotlib --user
pip install jupyter --user
```


Homebrew or Linuxbrew users
---------------------------

Please use [homebrew-ecell4](https://github.com/ecell/homebrew-ecell4)

https://github.com/ecell/homebrew-ecell4


Simple examples
---------------

Here are two extremely simple examples.
Please see http://ecell4.readthedocs.org for more details.

You need to add PYTHONPATH to import latest matplotlib only on Mac OSX.

### Running Python shell and E-Cell4

```shell
# on Windows or Linux
python

# on Mac
PYTHONPATH=~/Library/Python/2.7/lib/python/site-packages/ python
```

```
Python 2.7.6 (default, Mar 22 2014, 22:59:56)
[GCC 4.8.2] on linux2
Type "help", "copyright", "credits" or "license" for more information.
>>> from ecell4 import *
>>> sp = Species("B.A.C")
>>> print(sp.serial())
B.A.C
>>> print(unique_serial(sp))
A.B.C
```

### A reversible binding reaction

```python
%matplotlib inline
from ecell4 import *

with reaction_rules():
    A + B == C | (0.01, 0.3)

run_simulation(10, {'A': 60, 'B': 60}, solver='ode')
```

![png](https://raw.githubusercontent.com/ecell/ecell4/master/docs/output_7_0.png)
