================================
ecell4
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
