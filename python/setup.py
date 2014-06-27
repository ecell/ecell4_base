import sys
from distutils.core import setup, Command
from distutils.extension import Extension
import unittest

try:
    from Cython.Distutils import build_ext
except:
    print "You don't seem to have Cython installed. Please get a"
    print "copy from www.cython.org and install it"
    sys.exit(1)

class run_tests(Command):
    user_options = []

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        test_loader = unittest.defaultTestLoader
        suite = unittest.TestSuite()
        suite.addTest(test_loader.discover("tests/core"))
        # suite.addTest(test_loader.discover("tests/gillespie"))
        test_runner = unittest.TextTestRunner()
        test_runner.run(suite)

setup(
    packages = ["ecell4"],
    cmdclass = {'build_ext': build_ext, 'test': run_tests},
    ext_modules = [
        Extension("ecell4.core", sources=["ecell4/core.pyx"],
            include_dirs=["."], libraries=["ecell4-core"], language="c++"),
        Extension("ecell4.gillespie", sources=["ecell4/gillespie.pyx"],
            include_dirs=["."], libraries=["ecell4-core", "ecell4-gillespie"],
            language="c++")
        ])
