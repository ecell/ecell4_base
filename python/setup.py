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
        test_loader = unittest.TestLoader()
        def load_tests(module_name):
            return test_loader.discover(
                "tests/%s" % module_name, top_level_dir="tests")

        suite = unittest.TestSuite()
        suite.addTest(load_tests("core"))
        # suite.addTest(load_tests("gillespie"))
        # suite.addTest(load_tests("bd"))
        # suite.addTest(load_tests("ode"))
        # suite.addTest(load_tests("lattice"))
        suite.addTest(load_tests("reaction_reader"))
        # suite.addTest(load_tests("util"))
        test_runner = unittest.TextTestRunner()
        test_runner.run(suite)

setup(
    name = "ecell4",
    package_dir = {"": "lib"},
    package_data = {"ecell4.util": ["templates/*"]},
    packages = ["ecell4", "ecell4.util", "ecell4.reaction_reader"],
    cmdclass = {'build_ext': build_ext, 'test': run_tests},
    ext_modules = [
        Extension("ecell4.core", sources=["lib/ecell4/core.pyx"],
            include_dirs=["."], libraries=["ecell4-core"], language="c++"),
        Extension("ecell4.gillespie", sources=["lib/ecell4/gillespie.pyx"],
            include_dirs=["."], libraries=["ecell4-core", "ecell4-gillespie"],
            language="c++"),
        Extension("ecell4.bd", sources=["lib/ecell4/bd.pyx"],
            include_dirs=["."], libraries=["ecell4-core", "ecell4-bd"],
            language="c++"),
        Extension("ecell4.ode", sources=["lib/ecell4/ode.pyx"],
            include_dirs=["."], libraries=["ecell4-core", "ecell4-ode"],
            language="c++"),
        Extension("ecell4.lattice", sources=["lib/ecell4/lattice.pyx"],
            include_dirs=["."], libraries=["ecell4-core", "ecell4-lattice"],
            language="c++")
        ])
