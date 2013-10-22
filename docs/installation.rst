Installation
==================

Installing dependencies (Ubuntu Linux)
--------------------------------------------

::

   $ sudo aptitude install python
   $ sudo aptitude install pkg-config
   $ sudo aptitude install libgsl0-dev
   $ sudo aptitude install libboost1.53-all-dev
   $ sudo aptitude install libhdf5-serial-dev (please use Ubuntu newer than 12.10)
   $ sudo aptitude install python-pip
   $ pip install cython --user

Installing dependencies (Mac OSX)
-------------------------------------

::

   $ brew tap homebrew/science
   $ brew install pkg-config
   $ brew install gsl
   $ brew install boost
   $ brew install hdf5 --enable-cxx
   $ sudo easy_install pip
   $ pip install cython --user
   
Enviromental variables
--------------------------------

Please append these lines in .bashrc or .zshrc

in ubuntu

::

   export PATH=$HOME/.local/bin:$PATH
   export LD_LIBRARY_PATH=$HOME/ecell4/lib:$LD_LIBRARY_PATH


in MacOSX

::

   export PATH=$HOME/Library/Python/2.7/bin:$PATH
   export LD_LIBRARY_PATH=$HOME/ecell4/lib:$LD_LIBRARY_PATH


Installing E-Cell4 core
-----------------------------

core_python depends on core, you need to install core before you use core_python.

::

   $ PREFIX=$HOME/ecell4 bash install.sh core
   $ PREFIX=$HOME/ecell4 bash install.sh core_python

Installing E-Cell4 gillespie
----------------------------------

gillespie depends on core, you need to install core before you use gillespie.
gillespie_python depends on core_python, you need to install core_python before you use gillespie_python.

::

   $ PREFIX=$HOME/ecell4 bash install.sh gillespie
   $ PREFIX=$HOME/ecell4 PYTHONPATH=$PREFIX/lib/python2.7/site-packages bash install.sh gillespie_python

Installing E-Cell4 reaction_reader
---------------------------------------

reaction_reader depends on core_python, you need to install core_python before you use reaction_reader.

::

   $ PREFIX=$HOME/ecell4 bash install.sh reaction_reader


