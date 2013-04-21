from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

ext_modules = [
        Extension("PyEcell4",
            ["PyEcell4.pyx"],
            language = 'c++',
            libraries = ['ecell4-core'],
            )]

setup(
        name = "PyEcell4",
        cmdclass = {'build_ext' : build_ext},
        ext_modules = ext_modules
        )

#Gillespie
setup(
        name = "PyGillespie",
        cmdclass = {'build_ext' : build_ext},
        ext_modules = [
            Extension("PyGillespie",
                ["PyGillespie.pyx"],
                language = 'c++',
                libraries = ['ecell4-core', 'ecell4-gillespie'],
                )]
            )
