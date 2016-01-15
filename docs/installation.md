Installation
============

- [Docker users](#docker-users)
  - [Windows or Mac](#windows-or-mac)
  - [Linux](#linux)

- [Windows](#windows)
  - [Python2 series](#python2-series)
  - [Python3 series](#python3-series)
  
- [Mac](#mac)
  - [pip users](#pip-users)
  - [homebrew users](#homebrew-users)

- [Linux](#Linux)

- [Using E-Cell4 with jupyter](#using-e-cell4-with-jupyter)

- [Simple examples](#simple-examples)

Docker users
------------

If you have docker environment, you can try E-Cell4 easily.
You can pull E-Cell4 container with `docker pull ecell/ecell4`.

After the following steps, you should see Jupyter Notebook up and running (and E-Cell4 tutorials) in your web browser.

### Windows or Mac

1. Install [Docker Toolbox](https://www.docker.com/toolbox).
2. Run **Docker Quickstart Terminal**.
3. Run `docker run -d -p 443:8888 ecell/ecell4` in the terminal.
4. Open **192.168.99.100:443** with your favorite web browser.

### Linux

1. Install docker.
2. Run the following command.

    ```shell
    $ sudo docker pull ecell/ecell4
    $ sudo docker run -d -p 443:8888 ecell/ecell4
    ```

3. Open **localhost:443** with your favorite web browser.


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

### Python3 series

- [Python 3.5.1(**32bit**)](https://www.python.org/ftp/python/3.5.1/python-3.5.1.msi)
- [HDF5-1.8.16(**32-bit built with VS2015**)](http://www.hdfgroup.org/ftp/HDF5/current/bin/windows/extra/hdf5-1.8.16-win32-vs2015-shared.zip)

Please add python.exe, pip.exe path and `C:\Program Files (x86)\HDF_Group\HDF5\1.8.16\bin` to your **USER** PATH enviromental variable.
Next download numpy-1.10.4+vanilla-cp35-none-win32.whl and matplotlib-1.5.0-cp35-none-win32.whl from http://www.lfd.uci.edu/~gohlke/pythonlibs/
And run the following commands with command prompt.

```
pip install https://ci.appveyor.com/api/buildjobs/jpyueyasgwsannch/artifacts/python/dist/ecell4-4.0.0b2-cp35-none-win32.whl
pip install numpy-1.10.4+vanilla-cp35-none-win32.whl
pip install matplotlib-1.5.0-cp35-none-win32.whl
pip install -U jupyter
```

Mac
---

### pip users

1. Download [get-pip.py](https://bootstrap.pypa.io/get-pip.py)
2. Run the following commands
    ```shell
    sudo python get-pip.py
    # please select appropriate whl file for your Python version
    sudo pip install THEWHEELURL.whl
    # Mac default matplotlib is too old, you need to add these options to the pip command.
    pip install -U matplotlib --user
    sudo pip install -U jupyter
    ```

### homebrew users
Please see [homebrew-ecell4](https://github.com/ecell/homebrew-ecell4)

Linux
-----
Please use linuxbrew, see [homebrew-ecell4](https://github.com/ecell/homebrew-ecell4)

Using E-Cell4 with jupyter 
--------------------------

### Windows or Linux

```
jupyter-notebook
```

### Mac

You need to add user local Python site-package path to your PYTHONPATH to import latest matplotlib (instead of default matplotlib)
```
PYTHONPATH=~/Library/Python/2.7/lib/python/site-packages/ jupyter-notebook
```


Simple examples
---------------

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

![png](https://raw.githubusercontent.com/ecell/ecell4/master/docs/output_7_0.png)
