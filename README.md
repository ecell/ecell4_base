E-Cell System version 4
=======================

|   |  |
| ------------- | ------------- |
| Conda package for Windows, Mac and Linux | [![Conda Version](https://img.shields.io/conda/vn/conda-forge/ecell4_base.svg)](https://anaconda.org/conda-forge/ecell4_base) |
| PyPI (wheel) package for Linux | [![PyPI](https://img.shields.io/pypi/v/ecell4_base.svg)](https://pypi.python.org/pypi/ecell4_base) |

[![CircleCI](https://circleci.com/gh/ecell/ecell4_base.svg?style=svg)](https://circleci.com/gh/ecell/ecell4_base)
[![Documentation Status](https://readthedocs.org/projects/ecell4/badge/?version=latest)](http://ecell4.readthedocs.org/en/latest/?badge=latest)
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/87e076986e354b508f66af0a0ca3373d)](https://www.codacy.com/app/ecell/ecell4_base?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=ecell/ecell4_base&amp;utm_campaign=Badge_Grade)
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
Please refer to https://github.com/ecell/ecell4_docs

Quick start
-----------

Please refer to https://github.com/ecell/ecell4#quick-start

Installation
-------------

E-Cell4 does NOT support Python2.

ecell4_base does NOT support `pip install` on Windows and Mac. Please use `conda` instead.

### Windows, Mac, Linux

Install Miniconda with Python 3.7 for **64-bit** (from https://docs.conda.io/en/latest/miniconda.html)
and run these commands on your Terminal app.

```
conda install -c conda-forge ecell4_base
pip install ecell4
```

### Linux environment where you can NOT use conda (For example Google Colab)
We provide `ecell4_base` **wheel** package only for Linux.

```
python3 -m pip install ecell4 ecell4_base
```

### If you need to compile ecell4_base by yourself

Please refer to https://github.com/ecell/ecell4_base/blob/master/azure-pipelines.yml


Citation
========

If this package contributes to a project which leads to a scientific publication, I would appreciate a citation.

[![DOI](https://zenodo.org/badge/6348303.svg)](https://zenodo.org/badge/latestdoi/6348303)

Licensing terms
===============

This product is licensed under the terms of the [GNU General Public License v3](https://github.com/ecell/ecell4_base/blob/master/LICENSE),
See also [LICENSE](https://github.com/ecell/ecell4_base/blob/master/LICENSE) for the software included in this product.

- Copyright (c) 2010-, RIKEN

All rights reserved.
