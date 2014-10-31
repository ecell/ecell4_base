import sys
import glob
# from distutils.core import setup, Command
from setuptools import setup
from distutils.core import Command
from distutils.extension import Extension
import unittest

try:
    from Cython.Distutils import build_ext
except ImportError:
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
        suite.addTest(load_tests("util"))
        test_runner = unittest.TextTestRunner()
        test_runner.run(suite)

if sys.platform == "win32":
    dependent_libs = [
        'gsl', 'cblas', 'hdf5_cpp', 'hdf5']
    extra_compile_args = ["/EHsc", "/w"]
    with_cpp_shared_libraries = False
elif sys.platform == "darwin":
    dependent_libs = ['gsl', 'gslcblas', 'm', 'hdf5_cpp', 'hdf5']
    extra_compile_args = []
    # with_cpp_shared_libraries = False
    with_cpp_shared_libraries = False
else:
    dependent_libs = ['gsl', 'gslcblas', 'm', 'hdf5_cpp', 'hdf5']
    extra_compile_args = []
    # with_cpp_shared_libraries = False
    with_cpp_shared_libraries = True

if with_cpp_shared_libraries:
    ext_modules = [
        Extension("ecell4.core", sources=["lib/ecell4/core.pyx"],
            include_dirs=["."], libraries=["ecell4-core"], language="c++"),
        Extension("ecell4.egfrd", sources=["lib/ecell4/egfrd.pyx"],
            include_dirs=["."], libraries=["ecell4-core", "ecell4-egfrd"],
            language="c++", extra_compile_args=["-w"]),
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
            language="c++"),
        Extension("ecell4.meso", sources=["lib/ecell4/meso.pyx"],
            include_dirs=["."], libraries=["ecell4-core", "ecell4-meso"],
            language="c++"),
        ]
else:
    core_src = glob.glob("../ecell4/core/*.cpp")
    ext_modules = [
        Extension("ecell4.core", sources=["lib/ecell4/core.pyx"] + core_src,
            extra_compile_args=extra_compile_args,
            include_dirs=[".", ".."], libraries=dependent_libs, language="c++"),
        Extension("ecell4.egfrd",
            sources=["lib/ecell4/egfrd.pyx"]
                + glob.glob("../ecell4/egfrd/*.cpp") + core_src,
            extra_compile_args=extra_compile_args,
            libraries=dependent_libs, include_dirs=[".", ".."], language="c++"),
        Extension("ecell4.gillespie",
            sources=["lib/ecell4/gillespie.pyx"]
                + glob.glob("../ecell4/gillespie/*.cpp") + core_src,
            extra_compile_args=extra_compile_args,
            libraries=dependent_libs, include_dirs=[".", ".."], language="c++"),
        Extension("ecell4.bd",
            sources=["lib/ecell4/bd.pyx"]
                + glob.glob("../ecell4/bd/*.cpp") + core_src,
            extra_compile_args=extra_compile_args,
            libraries=dependent_libs, include_dirs=[".", ".."], language="c++"),
        Extension("ecell4.ode",
            sources=["lib/ecell4/ode.pyx"]
                + glob.glob("../ecell4/ode/*.cpp") + core_src,
            extra_compile_args=extra_compile_args,
            libraries=dependent_libs, include_dirs=[".", ".."], language="c++"),
        Extension("ecell4.lattice",
            sources=["lib/ecell4/lattice.pyx"]
                + glob.glob("../ecell4/lattice/*.cpp") + core_src,
            extra_compile_args=extra_compile_args,
            libraries=dependent_libs, include_dirs=[".", ".."], language="c++"),
        Extension("ecell4.meso",
            sources=["lib/ecell4/meso.pyx"]
                + glob.glob("../ecell4/meso/*.cpp") + core_src,
            extra_compile_args=extra_compile_args,
            libraries=dependent_libs, include_dirs=[".", ".."], language="c++"),
        ]

setup(
    name = "ecell4",
    package_dir = {"": "lib"},
    package_data = {"ecell4.util": [
        "templates/init_ipynb.js", "templates/init_cyjs.js", "templates/template.html",
        "templates/*.tmpl", "templates/ecelllogo/*.png"]},
    packages = ["ecell4",
        "ecell4.util", "ecell4.util.legacy"],
    cmdclass = {'build_ext': build_ext, 'test': run_tests},
    ext_modules = ext_modules
    )
