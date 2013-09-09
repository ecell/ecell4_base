================================
ecell4 |build-status|
================================

About
=====

E-Cell System, a software platform for modeling, simulation and analysis of complex, heterogeneous and multi-scale systems like the cell.

Dependencies
============

..

  $ sudo aptitude install libgsl0-dev libboost1.53-dev libhdf5-serial-dev


Install
=======

Do following instructions.
(To build Cython Wrapper it requires Cython 0.17 or later.)

..

  $ ./waf update --files="boost,doxygen" configure --prefix=${PREFIX} build install

Full installation (Read only)

..

  $ cd ${SRCPATH}

  for ode submodule:

  $ git clone git://github.com/headmyshoulder/odeint-v2 odeint-v2

  for spatiocyte submodule:

  $ cd ${SRCPATH}

  $ git clone git://github.com/ecell/ecell3.git ecell3

  $ cd ${SRCPATH}/ecell3

  $ ./autogen.sh && ./configure --prefix=${PREFIX} && make -j && make install

  $ cd ${SRCPATH}

  $ git clone git://github.com/ecell/ecell3-spatiocyte.git ecell3-spatiocyte

  $ cd ${SRCPATH}/ecell3-spatiocyte

  $ make

  install ecell4:

  $ cd ${SRCPATH}

  $ git clone git://github.com/ecell/ecell4.git ecell4

  $ cd ../ecell4

  $ PREFIX=/foo/bar \
  CPLUS_INCLUDE_PATH=${SRCPATH}/odeint-v2:${PREFIX}/include/ecell-3.2:\
  ${SRCPATH}/ecell3-spatiocyte:${SRCPATH}/epdp \
  LD_LIBRARY_PATH=${SRCPATH}/ecell3-spatiocyte:${SRCPATH}/epdp \
  LIBRARY_PATH=${SRCPATH}/ecell3-spatiocyte:${SRCPATH}/epdp \
  ./install.sh core bd gillespie ode egfrd spatiocyte

.. Build status badge
.. |build-status|
   image:: https://secure.travis-ci.org/ecell/ecell4.png
   :target: http://travis-ci.org/ecell/ecell4
   :alt: Build Status

Using homebrew (Mac OS X and OS X only)

The following commands will install core, gillespie, ode and reaction_reader

..

$ mv ecel4.rb /usr/local/Library/Formula

$ brew install ecell4

If you do not have pip or cython installed in your environment, homebrew will ask for your root password. If you do not want homebrew to do this, exit and do

for pip:

$ sudo easy_install pip

for cython:

$ sudo pip install cython

for cython without for local:

$ pip install cython --user

and add the install directory to your $PYTHONPATH

for csh and tcsh (assuming Python is version 2.7):

$ setenv PYTHONPATH $PYTHONPATH:/usr/local/lib/python2.7/site-packages

for other shells:

$ export PYTHONPATH=$PYTHONPATH:/usr/local/lib/python2.7/site-packages
