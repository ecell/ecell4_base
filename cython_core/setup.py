from distutils.core import setup
from Cython.Build import cythonize

'''
setup(
    name = "CoreObjects",
    ext_modules = cythonize('*.pyx'),
)
'''
setup(ext_modules = cythonize(
    "Species.pyx",
    language="c++",
    include_dirs="../core"
    ))
