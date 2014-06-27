import sys
from distutils.core import setup
from distutils.extension import Extension


try:
    from Cython.Distutils import build_ext
    except:
        print "You don't seem to have Cython installed. Please get a"
        print "copy from www.cython.org and install it"
        sys.exit(1)

setup(
    packages = ["ecell4"],
    cmdclass = {'build_ext': build_ext},
    ext_modules = [
        Extension("ecell4.core", sources=["ecell4/core.pyx"],
            include_dirs=["."], libraries=["ecell4-core"], language="c++"),
        Extension("ecell4.gillespie", sources=["ecell4/gillespie.pyx"],
            include_dirs=["."], libraries=["ecell4-core", "ecell4-gillespie"],
            language="c++")
        ])
