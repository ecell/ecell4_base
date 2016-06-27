import sys
import os.path
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

#XXX: $ cmake .
#XXX: $ make BesselTables
#XXX: $ cd python
#XXX: $ mkdir src
#XXX: $ cp ../ecell4 src/.
#XXX: $ cp ../ipynb src/.
#XXX: $ python setup.py sdict
#XXX: $ pip install dist/ecell-4.0.0.tar.gz
# STATIC_BUILD = True
STATIC_BUILD = False

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
with_egfrd = True
if "--disable-egfrd" in sys.argv:
    with_egfrd = False
    sys.argv.remove("--disable-egfrd")
    f = open('lib/ecell4/__init__.py', 'r+')
    contents = f.read()
    replaced_contents = contents.replace('from ecell4 import bd, ode, gillespie, egfrd, spatiocyte, meso', 'from ecell4 import bd, ode, gillespie, spatiocyte, meso')
    f.seek(0)
    f.write(replaced_contents)
    f.truncate()
    f.close()

if "--prefer-shared" in sys.argv:
    #XXX: This might be not a proper way to give a user defined parameter
    with_cpp_shared_libraries = True
    sys.argv.remove("--prefer-shared")

with_hdf5 = False
if "--hdf5" in sys.argv:
    #XXX: This might be not a proper way to give a user defined parameter
    with_hdf5 = True
    sys.argv.remove("--hdf5")

if sys.platform == "win32":
    with_hdf5 = True  #XXX: forced
    if sys.version_info.major == 2:
        dependent_libs = ['gsl', 'cblas']
        extra_compile_args = ["/EHsc", "-DHAVE_CONFIG_H", "-DHAVE_INLINE"]
    elif sys.version_info.major == 3:
        dependent_libs = ['gsl', 'gslcblas']
        extra_compile_args = ["/EHsc", "/w", "-DHAVE_CONFIG_H"]  # "-DHAVE_INLINE"
    # extra_compile_args.append('-DNO_BESSEL_TABLE')
elif sys.platform == "darwin":
    with_hdf5 = True  #XXX: forced
    dependent_libs = ['gsl', 'gslcblas', 'm']
    extra_compile_args = ["-DNO_BESSEL_TABLE", "-DHAVE_CONFIG_H"]
else: # for linux
    if not with_cpp_shared_libraries:
        with_hdf5 = True  #XXX: forced
    dependent_libs = ['gsl', 'gslcblas', 'm']
    extra_compile_args = ["-DNO_BESSEL_TABLE", "-DHAVE_CONFIG_H"]

if "--disable-hdf5" in sys.argv:
    with_hdf5 = False
    sys.argv.remove("--disable-hdf5")

if with_hdf5:
    dependent_libs.extend(['hdf5_cpp', 'hdf5'])
    extra_compile_args.append("-DWITH_HDF5")
    if sys.platform == "win32":
        extra_compile_args.extend(
            ["-D_HDF5USEDLL_", "-DHDF5CPP_USEDLL", "-DH5_BUILT_AS_DYNAMIC_LIB"])

src_path = ".."

if STATIC_BUILD:
    src_path = "./src"
    with_cpp_shared_libraries = False
    if not (os.path.isdir(os.path.join(src_path, "ecell4"))
            and os.path.isdir(os.path.join(src_path, "ipynb"))):
        raise IOError(
            ("The source directory '{0}/ecell4' or '{0}/ipynb' was not found."
             + "Do as follows: cp -r ../ecell4 {0}/.").format(src_path))
    elif (os.path.islink(os.path.join(src_path, "ecell4"))
            or os.path.islink(os.path.join(src_path, "ipynb"))):
        raise UserWarning(
            ("The source directory '{0}/ecell4' or '{0}/ipynb' is a symbolic link."
             + "It might fail to be scanned properly."))

if with_cpp_shared_libraries:
    ext_modules = [
        Extension("ecell4.core", sources=["lib/ecell4/core.pyx"],
            include_dirs=["."], libraries=["ecell4-core"], language="c++",
            extra_compile_args=extra_compile_args),
        Extension("ecell4.egfrd", sources=["lib/ecell4/egfrd.pyx"],
            include_dirs=["."], libraries=["ecell4-core", "ecell4-egfrd"],
            language="c++", extra_compile_args=extra_compile_args + ["-w"]),
        Extension("ecell4.gillespie", sources=["lib/ecell4/gillespie.pyx"],
            include_dirs=["."], libraries=["ecell4-core", "ecell4-gillespie"],
            language="c++", extra_compile_args=extra_compile_args),
        Extension("ecell4.bd", sources=["lib/ecell4/bd.pyx"],
            include_dirs=["."], libraries=["ecell4-core", "ecell4-bd"],
            language="c++", extra_compile_args=extra_compile_args),
        Extension("ecell4.ode", sources=["lib/ecell4/ode.pyx"],
            include_dirs=["."], libraries=["ecell4-core", "ecell4-ode"],
            language="c++", extra_compile_args=extra_compile_args),
        Extension("ecell4.spatiocyte", sources=["lib/ecell4/spatiocyte.pyx"],
            include_dirs=["."], libraries=["ecell4-core", "ecell4-spatiocyte"],
            language="c++", extra_compile_args=extra_compile_args),
        Extension("ecell4.meso", sources=["lib/ecell4/meso.pyx"],
            include_dirs=["."], libraries=["ecell4-core", "ecell4-meso"],
            language="c++", extra_compile_args=extra_compile_args),
        ]
    ext_modules = cythonize(ext_modules)
else:
    core_src = glob.glob(os.path.join(src_path, "ecell4/core/*.cpp"))

    ext_modules = [
        Extension("ecell4.core", sources=["lib/ecell4/core.pyx"] + core_src,
            extra_compile_args=extra_compile_args,
            include_dirs=[src_path], libraries=dependent_libs, language="c++"),
        Extension("ecell4.gillespie",
            sources=["lib/ecell4/gillespie.pyx"]
                + glob.glob(os.path.join(src_path, "ecell4/gillespie/*.cpp")) + core_src,
            extra_compile_args=extra_compile_args,
            libraries=dependent_libs, include_dirs=[src_path], language="c++"),
        Extension("ecell4.bd",
            sources=["lib/ecell4/bd.pyx"]
                + glob.glob(os.path.join(src_path, "ecell4/bd/*.cpp")) + core_src,
            extra_compile_args=extra_compile_args,
            libraries=dependent_libs, include_dirs=[src_path], language="c++"),
        Extension("ecell4.ode",
            sources=["lib/ecell4/ode.pyx"]
                + glob.glob(os.path.join(src_path, "ecell4/ode/*.cpp")) + core_src,
            extra_compile_args=extra_compile_args,
            libraries=dependent_libs, include_dirs=[src_path], language="c++"),
        Extension("ecell4.spatiocyte",
            sources=["lib/ecell4/spatiocyte.pyx"]
                + glob.glob(os.path.join(src_path, "ecell4/spatiocyte/*.cpp")) + core_src,
            extra_compile_args=extra_compile_args,
            libraries=dependent_libs, include_dirs=[src_path], language="c++"),
        Extension("ecell4.meso",
            sources=["lib/ecell4/meso.pyx"]
                + glob.glob(os.path.join(src_path, "ecell4/meso/*.cpp")) + core_src,
            extra_compile_args=extra_compile_args,
            libraries=dependent_libs, include_dirs=[src_path], language="c++"),
        ]
    if with_egfrd:
        ext_modules.append(Extension("ecell4.egfrd",
            sources=["lib/ecell4/egfrd.pyx"]
                + glob.glob(os.path.join(src_path, "ecell4/egfrd/*.cpp")) + core_src,
            extra_compile_args=extra_compile_args,
            libraries=dependent_libs, include_dirs=[src_path], language="c++"))

    ext_modules = cythonize(ext_modules)

setup(
    name = "ecell",
    version = "4.0.1",
    package_dir = {"": "lib"},
    package_data = {"ecell4.util": [
        "templates/init_ipynb.js", "templates/init_cyjs.js", "templates/template.html",
        "templates/*.tmpl", "templates/ecelllogo/*.png"]},
    data_files = [('ecell4ipynb/Licenses', glob.glob(os.path.join(src_path, 'licenses/*')))],
    # data_files = [('ecell4ipynb/Licenses', glob.glob(os.path.join(src_path, 'licenses/*'))),
    #               ('ecell4ipynb', [os.path.join(src_path, 'ipynb/index.ipynb')]),
    #               ('ecell4ipynb/Tutorials', glob.glob(os.path.join(src_path, 'ipynb/Tutorials/*.ipynb'))),
    #               ('ecell4ipynb/Examples', glob.glob(os.path.join(src_path, 'ipynb/Examples/*.ipynb'))),
    #               ('ecell4ipynb/Tests', glob.glob(os.path.join(src_path, 'ipynb/Tests/*.ipynb'))),
    #               ('ecell4ipynb/Sandbox', glob.glob(os.path.join(src_path, 'ipynb/Sandbox/*.ipynb'))),
    #               ],
    packages = ["ecell4", "ecell4.util", "ecell4.extra"],
    cmdclass = {'build_ext': build_ext, 'test': run_tests},
    license = "the GNU General Public License v2",
    author = "Kazunari Kaizu",
    author_email = "kaizu@riken.jp",
    url = "https://github.com/ecell/ecell4",
    ext_modules = ext_modules
    )
