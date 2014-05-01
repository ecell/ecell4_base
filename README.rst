================================
E-Cell System version 4 
================================

What is E-Cell System?
=======================

E-Cell System, a software platform for modeling, simulation and analysis of complex, heterogeneous and multi-scale systems like the cell.

Requirements
=============

..

  $ sudo apt-get install libgsl0-dev libboost-dev libboost-test-dev libboost-regex-dev libhdf5-serial-dev

  $ sudo apt-get instal python-dev cython

How to install?
================

Do following instructions on Ubuntu 14.04 LTS (Trusty Tahr).

..

   $ wget https://github.com/ecell/ecell4/archive/master.zip
   
   $ unzip master.zip
   
   $ cd ecell4-master
   
   $ PREFIX=/path/to PYTHONPATH=/path/to/lib/python2.7/site-packages ./install.sh

How to use?
============

..

   $ LD_LIBRARY_PATH=/pat/to/lib PYTHONPATH=/path/to/lib/python2.7/site-packages python
   Python 2.7.6 (default, Mar 22 2014, 22:59:56) 
   [GCC 4.8.2] on linux2
   Type "help", "copyright", "credits" or "license" for more information.
   >>> from ecell4.core import *
   >>> sp = Species("B.A.C")
   >>> print sp.serial()
   A.B.C
   >>> 

.. Build status badge
.. |build-status|
   image:: https://secure.travis-ci.org/ecell/ecell4.png
   :target: http://travis-ci.org/ecell/ecell4
   :alt: Build Status
