E-Cell System version 4
=======================

| Theme | Status |
| ------------- | ------------- |
| Conda package for Windows, Mac and Linux | [![Conda Version](https://img.shields.io/conda/vn/conda-forge/ecell4_base.svg)](https://github.com/conda-forge/ecell4_base-feedstock) |
| PyPI (wheel) package for Linux | [![PyPI](https://img.shields.io/pypi/v/ecell4_base.svg)](https://pypi.python.org/pypi/ecell4_base) |

[![Documentation Status](https://readthedocs.org/projects/ecell4/badge/?version=latest)](http://ecell4.readthedocs.org/en/latest/?badge=latest)

What is E-Cell System?
----------------------

E-Cell System is a software platform for modeling, simulation and analysis of complex, heterogeneous and multi-scale systems like a cell.

See https://ecell4.e-cell.org/

Build macOS arm64 wheel
------------------------

The following steps will build:
`dist/ecell4_base-2.1.2-cp310-cp310-macosx_26_0_arm64.whl`
on macOS Tahoe (Apple Silicon, arm64).

```bash
# 1) Install required build dependencies
brew install cmake hdf5 boost gsl

# 2) Build wheel for arm64
cd /path/to/ecell4_base
rm -rf dist build
uv venv --python 3.10
# or
# uv venv --python 3.11
# uv venv --python 3.12
# uv venv --python 3.13
# uv venv --python 3.14
source .venv/bin/activate
uv pip install pip
CMAKE_OSX_ARCHITECTURES=arm64 pip wheel . -w dist --no-deps

# 3) (Optional) Verify install from the built wheel in the above venv
pip install dist/ecell4_base-2.1.2-cp310-cp310-macosx_26_0_arm64.whl
python -c "from ecell4_base.core import *; sp1 = Species('A', 0.0025, 1); print(sp1.list_attributes())"
```
