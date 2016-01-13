import sys
import glob
import unittest

from setuptools import setup
from distutils.core import Command, Extension
# from distutils.core import setup, Command, Extension

try:
    from Cython.Build import cythonize
    from Cython.Distutils import build_ext
except ImportError:
    print("You don't seem to have Cython installed. Please get a")
    print("copy from www.cython.org and install it")
    sys.exit(1)

sys.path.append("./lib")

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
        # suite.addTest(load_tests("spatiocyte"))
        suite.addTest(load_tests("util"))
        test_runner = unittest.TextTestRunner()
        test_runner.run(suite)

with_cpp_shared_libraries = False
if "--prefer-shared" in sys.argv:
    #XXX: This might be not a proper way to give a user defined parameter
    with_cpp_shared_libraries = True
    sys.argv.remove("--prefer-shared")

if sys.platform == "win32":
    if sys.version_info.major == 2:
        dependent_libs = ['gsl', 'cblas', 'hdf5_cpp', 'hdf5']
        extra_compile_args = ["/EHsc", "/w", "-DHAVE_CONFIG_H", "-DHAVE_INLINE", "-DWITH_HDF5", "-D_HDF5USEDLL_", "-DHDF5CPP_USEDLL", "-DH5_BUILT_AS_DYNAMIC_LIB"]
    elif sys.version_info.major == 3:
        # dependent_libs = ['gsl', 'gslcblas', 'hdf5_cpp', 'hdf5']
        #extra_compile_args = ["/EHsc", "/w", "-DHAVE_CONFIG_H", "-DWITH_HDF5", "-D_HDF5USEDLL_", "-DHDF5CPP_USEDLL", "-DH5_BUILT_AS_DYNAMIC_LIB"]  # "-DHAVE_INLINE"
        dependent_libs = ['gsl', 'gslcblas']
        extra_compile_args = ["/EHsc", "/w", "-DHAVE_CONFIG_H"]
elif sys.platform == "darwin":
    dependent_libs = ['gsl', 'gslcblas', 'm', 'hdf5_cpp', 'hdf5']
    extra_compile_args = ["-DWITH_HDF5"]
else: # for linux
    dependent_libs = ['gsl', 'gslcblas', 'm', 'hdf5_cpp', 'hdf5']
    extra_compile_args = []

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
        Extension("ecell4.spatiocyte", sources=["lib/ecell4/spatiocyte.pyx"],
            include_dirs=["."], libraries=["ecell4-core", "ecell4-spatiocyte"],
            language="c++"),
        Extension("ecell4.meso", sources=["lib/ecell4/meso.pyx"],
            include_dirs=["."], libraries=["ecell4-core", "ecell4-meso"],
            language="c++"),
        ]
    ext_modules = cythonize(ext_modules)
else:
    # import subprocess
    # import os.path
    # path_to_egfrd = os.path.join(os.path.abspath(".."), "ecell4", "egfrd")
    # sjy_table_path = os.path.join(path_to_egfrd, "SphericalBesselTable.hpp")
    # cjy_table_path = os.path.join(path_to_egfrd, "CylindricalBesselTable.hpp")

    # ext_modules = cythonize(
    #         glob.glob("lib/ecell4/*.pyx"),
    #         sources=glob.glob("../ecell4/*/*.cpp"),  #XXX: Not safe
    #         include_path=[".", ".."],
    #         extra_compile_args=extra_compile_args,
    #         libraries=dependent_libs, language="c++", setup_requires=["Cython"])

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
        Extension("ecell4.spatiocyte",
            sources=["lib/ecell4/spatiocyte.pyx"]
                + glob.glob("../ecell4/spatiocyte/*.cpp") + core_src,
            extra_compile_args=extra_compile_args,
            libraries=dependent_libs, include_dirs=[".", ".."], language="c++"),
        Extension("ecell4.meso",
            sources=["lib/ecell4/meso.pyx"]
                + glob.glob("../ecell4/meso/*.cpp") + core_src,
            extra_compile_args=extra_compile_args,
            libraries=dependent_libs, include_dirs=[".", ".."], language="c++"),
        ]
    ext_modules = cythonize(ext_modules)

setup(
    name = "ecell4",
    version = "4.0.0b2",
    package_dir = {"": "lib"},
    package_data = {"ecell4.util": [
        "templates/init_ipynb.js", "templates/init_cyjs.js", "templates/template.html",
        "templates/*.tmpl", "templates/ecelllogo/*.png"]},
    data_files = [('ecell4ipynb', ['../ipynb/index.ipynb']),
                  ('ecell4ipynb/Tutorials', glob.glob('../ipynb/Tutorials/*.ipynb'))],
    packages = ["ecell4", "ecell4.util", "ecell4.extra"],
    cmdclass = {'build_ext': build_ext, 'test': run_tests},
    ext_modules = ext_modules
    )
