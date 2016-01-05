Installation
============

- [Docker users](#docker-users)
  - [Windows or Mac](#windows-or-mac)
  - [Linux](#linux)

- [Windows](#windows)
- [Mac or Linux](#mac-or-linux)

- [Using E-Cell4 with jupyter](#using-e-cell4-with-jupyter)

- [Simple examples](#simple-examples)

Docker users
------------

If you have docker environment, you can try E-Cell4 easily.
You can pull E-Cell4 container with `docker pull ecell/ecell4`.

### Windows or Mac

1. Install [Docker Toolbox](https://www.docker.com/toolbox).
2. Run **Docker Quickstart Terminal**.
3. Run `docker run -d -p 443:8888 ecell/ecell4` in the terminal.
4. Open **192.168.99.100:443** with your favorite web browser.  
5. You should see Jupyter Notebook up and running (and E-Cell4 tutorials) in your web browser.

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

- [Python 2.7.11(**32bit**)](https://www.python.org/ftp/python/2.7.11/python-2.7.11.msi)
- [HDF5-1.8.16 Pre-built Binary(**32-bit**)](http://www.hdfgroup.org/ftp/HDF5/current/bin/windows/extra/hdf5-1.8.16-win32-vs2012-shared.zip)
- [Visual C++ Compiler for Python 2.7](http://aka.ms/vcpython27)

Please add `C:\Python27`, `C:\Python27\Scripts` and `C:\Program Files (x86)\HDF_Group\HDF5\1.8.16\bin` to your **PATH** enviromental variable.

And run the following commands with command prompt.
matplotlib depends on numpy. It takes some time to build numpy, please be patient.
Although jupyter is optional, we strongly recommend that you run E-Cell4 with jupyter.

```shell
pip install https://ci.appveyor.com/api/buildjobs/aju5rykh88bb88ns/artifacts/python/dist/ecell4-4.0.0b2-cp27-none-win32.whl
pip install matplotlib
pip install jupyter
```


Mac or Linux
------------

Please use [homebrew-ecell4](https://github.com/ecell/homebrew-ecell4)

https://github.com/ecell/homebrew-ecell4

### Note about Jupyter and matplotlib for Mac

```shell
sudo python get-pip.py
# Mac default matplotlib is too old, so you need to add these options to the pip command.
# These commands can be used not only for Mac but also for Linux.
pip install -U matplotlib --user
sudo pip install jupyter
```

Using E-Cell4 with jupyter 
--------------------------

### Windows or Linux

```
jupyter-notebook
```

### Mac

You need to add PYTHONPATH to import latest matplotlib (instead of default matplotlib)
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
