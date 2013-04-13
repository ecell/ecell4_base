from distutils.core import setup
from Cython.Build import cythonize

# core
setup(ext_modules = cythonize(
    "PySpecies.pyx",
    language="c++",
    include_dirs="../core"
    ))

setup(ext_modules = cythonize(
    "PyReactionRule.pyx",
    language="c++",
    include_dirs="../core"
    ))

setup(ext_modules = cythonize(
    "PyCompartmentSpace.pyx",
    language="c++",
    include_dirs="../core"
    ))

setup(ext_modules = cythonize(
    "PyNetworkModel.pyx",
    language="c++",
    include_dirs="../core"
    ))

setup(ext_modules = cythonize(
    "PyODEWorld.pyx",
    language="c++",
    include_dirs="../ode"
    ))
