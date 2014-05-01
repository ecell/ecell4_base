================================
ecell4 
================================

About
=====

E-Cell System, a software platform for modeling, simulation and analysis of complex, heterogeneous and multi-scale systems like the cell.

Dependencies
============

..

  $ sudo apt-get install libgsl0-dev libboost-dev libboost-test-dev libboost-regex-dev libhdf5-serial-dev

  $ sudo apt-get instal python-dev cython

Install
=======

Do following instructions on Ubuntu 14.04 LTS (Trusty Tahr).

..

   $ wget https://github.com/ecell/ecell4/archive/master.zip
   
   $ unzip master.zip
   
   $ cd ecell4-master
   
   $ PREFIX=/path/to PYTHONPATH=/path/to/lib/python2.7/site-packages ./install.sh

.. Build status badge
.. |build-status|
   image:: https://secure.travis-ci.org/ecell/ecell4.png
   :target: http://travis-ci.org/ecell/ecell4
   :alt: Build Status
