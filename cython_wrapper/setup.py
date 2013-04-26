from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

ext_modules = [
        Extension("Ecell4",
            ["Ecell4.pyx"],
            language = 'c++',
            libraries = ['ecell4-core'],
            )]

setup(
        name = "Ecell4",
        cmdclass = {'build_ext' : build_ext},
        ext_modules = ext_modules
        )

#Gillespie
setup(
        name = "Gillespie",
        cmdclass = {'build_ext' : build_ext},
        ext_modules = [
            Extension("Gillespie",
                ["Gillespie.pyx"],
                language = 'c++',
                libraries = ['ecell4-core', 'ecell4-gillespie', 'hdf5', 'hdf5_cpp'],
                )]
            )
