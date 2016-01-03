E-Cell System version 4
=======================

[![Build Status](https://travis-ci.org/ecell/ecell4.svg?branch=develop)](https://travis-ci.org/ecell/ecell4)
[![Gitter](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/ecell/ecell4?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge)

## What is E-Cell System?

E-Cell System, a software platform for modeling, simulation and analysis of complex, heterogeneous and multi-scale systems like the cell.

- [Dockerized E-Cell4 Jupyter notebooks](#dockerized-e-cell4-jupyter-notebooks)
    - [For Windows or Mac](#for-windows-or-mac)
    - [For Linux](#for-linux)

- [Native binary installation](#native-binary-installation)
    - [Windows](#windows)
    - [Mac or Linux](#mac-or-linux)

- [Simple examples](#simple-examples)

Dockerized E-Cell4 Jupyter notebooks
------------------------------------

If you have docker environment, you can try E-Cell4 easily.
You can pull E-Cell4 container with `docker pull ecell/ecell4`

### For Windows or Mac

1. Install [Docker Toolbox](https://www.docker.com/toolbox)
2. Run Kitematic
3. Search container with **ecell4**, and create ecell4 container

  ![png](https://raw.githubusercontent.com/ecell/ecell4/master/docs/kitematic1.png)

4. Click **WEB PREVIEW** button in **HOME** tab

  ![png](https://raw.githubusercontent.com/ecell/ecell4/develop/docs/images/kitematic2.png)
  
5. You should see Jupyter Notebook up and running (and E-Cell4 tutorials) in your default browser.

### For Linux

```shell
$ sudo docker pull ecell/ecell4
$ sudo docker run -d -p 443:8888 ecell/ecell4
```

You'll be able to E-Cell4 notebooks at http://THE_IP_RUNNING_DOCKER:443


Native binary installation
--------------------------

### Windows

#### Requirements

Please use 32bit Python, even if you use 64bit Windows.
We have NOT supported 64bit Python yet.

- [Python 2.7.11(**32bit**)](https://www.python.org/ftp/python/2.7.11/python-2.7.11.msi)
- HDF5-1.8.16 Pre-built Binary(**32-bit**) http://www.hdfgroup.org/ftp/HDF5/current/bin/windows/extra/hdf5-1.8.16-win32-vs2015-shared.zip

Please add `C:\Python27`, `C:\Python27\Scripts` and `C:\Program Files (x86)\HDF_Group\HDF5\1.8.16\bin` to your **PATH** enviromental variable.

And run following command with command prompt.
```
pip install https://ci.appveyor.com/api/buildjobs/aju5rykh88bb88ns/artifacts/python/dist/ecell4-4.0.0b2-cp27-none-win32.whl
```

#### Jupyter for Windows
We recommend you run E-Cell4 models from Jupyter notebook.
Below is Jupyter notebook(and matplotlib) installation for Windows.

- Install [Visual C++ Compiler for Python 2.7](http://aka.ms/vcpython27)
- Install Jupyter notebook and matplotlib

  ```
  pip install -U jupyter
  pip install matplotlib
  ```

matplotlib depends on numpy. It takes some time to build numpy, please be patient.

### Mac or Linux

Please use [homebrew-ecell4](https://github.com/ecell/homebrew-ecell4)

https://github.com/ecell/homebrew-ecell4

#### Jupyter for Mac or Linux

We recommend you run E-Cell4 models from Jupyter notebook.
Below is Jupyter notebook(and matplotlib) installation for Mac or Linux.

```shell
sudo python get-pip.py
sudo pip install matplotlib
sudo pip install jupyter
```

Simple examples
---------------

Here are two extremely simple examples.
Please see http://ecell4.readthedocs.org for more details.

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

#### A reversible binding reaction

```python
%matplotlib inline
from ecell4 import *

with reaction_rules():
    A + B == C | (0.01, 0.3)

run_simulation(10, {'A': 60, 'B': 60}, solver='ode')
```

![png](https://raw.githubusercontent.com/ecell/ecell4/master/docs/output_7_0.png)


