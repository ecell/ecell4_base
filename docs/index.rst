.. raw:: html

   <img src="https://raw.githubusercontent.com/ecell/ecell4/master/docs/images/ecell-logo-with-title.png">

.. raw:: html


Welcome to E-Cell System version 4 documentation!
=================================================

What is E-Cell System?
----------------------

E-Cell System is a software platform for modeling, simulation and
analysis of complex, heterogeneous and multi-scale systems like the
cell.

E-Cell4 is a free and open-source software licensed under the GNU
General Public License version 2.

The source code is available on
`GitHub <https://github.com/ecell/ecell4>`__.

For installation instructions, please see :doc:`installation`. E-Cell4 runs on Windows, Mac OSX and Linux.

4.0.1 (Released May 20, 2016)
-----------------------------------

Enhancement
~~~~~~~~~~~~~~

- Progressbar on Jupyter Notebook
- Better styles for visualization especially for matplotlib
- A member function to distinguish Species for molecules and structures
- Enable user-defined macros for Python decorators
- Enable to save a movie as mp4
- Enable to set tolerances for ODE through a Factory class
- 2D visualization of 3D world for animation
- getVelocity for ODE simulation
- Implicit dissolution of bindings in the rule-based modeling
- Higher resolution for plot_movie
- Enable to save time course data of NumberObservers
- FixedIntervalTrajectoryObserver should have its small interval to check the periodic boundary effect separately from logging interval
- Convenient functions for declaring fundamental objects
- Enable to set ranges to show trajectories on viz.plot_trajectory

Bug fix
~~~~~~~~~

- Fixed kernel death that would occur when there is no species in the viz.plot_movie_with_matplotlib(..., species_list)
- Implement egfrd::World::list_species()
- Spatiocyte doesn't accept the second-order reaction with translocation
- Wrong absorption rate for mesoscopic simulations


Getting Started
---------------

.. toctree::
   :glob:

   installation

Tutorials
---------

We assume that you run the following tutorials from Jupyter Notebook.
Please note that magic functions like ``%matplotlib`` work only on Jupyter Notebook.

.. toctree::
   :glob:
   
   tutorial1/tutorial1
   tutorial2/tutorial2
   tutorial3/tutorial3
   tutorial4/tutorial4
   tutorial5/tutorial5
   tutorial6/tutorial6
   tutorial7/tutorial7
   tutorial8/tutorial8
   tutorial9/tutorial9
   tutorial10/tutorial10

API
---

.. toctree::
   :glob:

   api/core
   api/gillespie
   api/ode
   api/meso
   api/spatiocyte
   api/bd
   api/egfrd
   api/util
   api/util.viz
   api/util.decorator


Examples
--------

.. toctree::
   :glob:
   
   Attractors/Attractors
   Drosophila/Drosophila
   DualPhospho/DualPhospho
   Glycolysis/Glycolysis
   HodgkinHuxley/HodgkinHuxley
   SimpleEquilibrium/SimpleEquilibrium
   MinDEMeso/MinDEMeso
   MinDESpatiocyte/MinDESpatiocyte
   
