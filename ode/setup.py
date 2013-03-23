from distutils.core import setup
from Cython.Build import cythonize

setup(
    name = "ODEApp",
    ext_modules = cythonize('*.pyx'),
)
