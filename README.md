E-Cell System version 4
=======================

[![Build Status](https://travis-ci.org/ecell/ecell4-base.svg?branch=develop)](https://travis-ci.org/ecell/ecell4-base)
[![CircleCI](https://circleci.com/gh/ecell/ecell4-base.svg?style=svg)](https://circleci.com/gh/ecell/ecell4-base)
[![Build status](https://ci.appveyor.com/api/projects/status/github/ecell/ecell4?svg=true)](https://ci.appveyor.com/project/kaizu/ecell4)
[![Documentation Status](https://readthedocs.org/projects/ecell4/badge/?version=latest)](http://ecell4.readthedocs.org/en/latest/?badge=latest)
[![PyPI](https://img.shields.io/pypi/v/ecell.svg)](https://pypi.python.org/pypi/ecell)
[![License: GPL v2](https://img.shields.io/badge/license-GPL%20v2-blue.svg)](https://github.com/ecell/ecell4-base/blob/master/licenses/LICENSE)
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/87e076986e354b508f66af0a0ca3373d)](https://www.codacy.com/app/ecell/ecell4-base?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=ecell/ecell4-base&amp;utm_campaign=Badge_Grade)
[![Slack Status](https://img.shields.io/badge/chat-on%20slack-50baa6.svg)](https://ecell-project.herokuapp.com/)
<!---[![Slack Status](https://ecell-project.herokuapp.com/badge.svg)](https://ecell-project.herokuapp.com/)--->

What is E-Cell System?
----------------------

E-Cell System is a software platform for modeling, simulation and analysis of complex, heterogeneous and multi-scale systems like a cell.

E-Cell4 accepts multi-algorithms, multi-timescales and multi-spatial-representations as its central feature.

Features
--------

- Single particle simulations, i.e. [The enhanced Green's Function Reaction Dynamics (eGFRD) method](http://gfrd.org), [Spatiocyte](http://spatiocyte.org) (a lattice-based method), and the Reaction Brownian Dynamics (RBD) method
- Ordinary differential equations, Gillespie algorithm (the direct method), and spatial Gillespie algorithm (the next subvolume method)
- Rule-based modeling
- Python programmable

Try online
----------

You can try this package online with Google Colaboratory.
Please refer to the https://github.com/ecell/ecell4-docs

Installation
-------------

Please see [INSTALL.md](https://github.com/ecell/ecell4-base/blob/master/INSTALL.md).
Basically you can install E-Cell4 on any OS just by running
```
pip install ecell4
```

Simple examples
---------------

Here are two extremely simple examples, See http://ecell4.readthedocs.org for more details on running E-Cell4.

```
Python 2.7.6 (default, Mar 22 2014, 22:59:56)
[GCC 4.8.2] on linux2
Type "help", "copyright", "credits" or "license" for more information.
>>> from ecell4.core import *
>>> sp = Species("A.B.C")
>>> print sp.serial()
A.B.C
>>>
```

### Binding and unbinding reactions

```python
%matplotlib inline
from ecell4 import *

with reaction_rules():
    A + B == C | (0.01, 0.3)

run_simulation(10, {'A': 60, 'B': 60})
```

![png](./readthedocs/images/output_7_0.png)

### Diffusion on a spherical surface

```python
%matplotlib inline
from ecell4 import *

with species_attributes():
    A | {'D': '1', 'location': 'M'}

surface = Sphere(ones() * 0.5, 0.5).surface()
obs = FixedIntervalTrajectoryObserver(1e-4)
run_simulation(
    0.4, y0={'A': 10}, structures={'M': surface},
    solver='spatiocyte', observers=obs, return_type=None)

viz.plot_trajectory(obs, interactive=False)
```

![png](./readthedocs/images/hairball.png)

Tutorials, Examples, and API (Static Site)
==========================================
https://ecell4.readthedocs.io/en/latest/

Citation
========

If this package contributes to a project which leads to a scientific publication, I would appreciate a citation.

[![DOI](https://zenodo.org/badge/6348303.svg)](https://zenodo.org/badge/latestdoi/6348303)

Licensing terms
===============

This product is licensed under the terms of the [GNU General Public License v2](https://github.com/ecell/ecell4/blob/master/licenses/LICENSE),
See [NOTICE](https://github.com/ecell/ecell4/blob/master/licenses/NOTICE.txt) for the software included in this product.

- Copyright (c) 2010-, RIKEN

All rights reserved.
