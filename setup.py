import os
import re
import sys
import platform
import sysconfig
import subprocess

from setuptools import setup, Extension
from setuptools.command.install import install
from distutils.core import Command
from distutils.version import LooseVersion
import shutil

from logging import getLogger, StreamHandler, Formatter, DEBUG
formatter = Formatter('%(levelname)s:%(name)s: %(message)s')
handler = StreamHandler()
handler.setLevel(DEBUG)
handler.setFormatter(formatter)
logger = getLogger(__name__)
logger.setLevel(DEBUG)
logger.addHandler(handler)

try:
    from Cython.Build import cythonize
    from Cython.Distutils import build_ext
except ImportError:
    logger.error("You don't seem to have Cython installed. Please get a")
    logger.error("copy from www.cython.org and install it")
    sys.exit(1)


class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=''):
        super(CMakeExtension, self).__init__(name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)


shared_libs_dir = os.path.abspath(os.path.join('build', 'shared_libs'))

class CMakeBuild(build_ext):
    def build_extensions(self):
        try:
            out = subprocess.check_output(['cmake', '--version'])
        except OSError:
            raise RuntimeError("CMake must be installed to build the following extensions: " +
                               ", ".join(e.name for e in self.extensions))

        cmake_minimum_version_for_windows = '3.1.0'
        if platform.system() == "Windows":
            cmake_version = LooseVersion(re.search(r'version\s*([\d.]+)', out.decode()).group(1))
            if cmake_version < cmake_minimum_version_for_windows:
                raise RuntimeError("CMake >= {} is required on Windows".format(
                    cmake_minimum_version_for_windows))

        for ext in self.extensions:
            if isinstance(ext, CMakeExtension):
                self.build_cmake_extension(ext)
            else:
                self.build_extension(ext)

    def build_extension(self, ext):
        ext.sources = list(map(os.path.relpath, ext.sources))
        ext.include_dirs += [
            self.build_temp
        ]
        ext.library_dirs = [shared_libs_dir]
        super(CMakeBuild, self).build_extension(ext)

    def build_cmake_extension(self, ext):
        # extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))
        extdir = shared_libs_dir
        cmake_args = ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=' + extdir]
                      # '-DPYTHON_EXECUTABLE=' + sys.executable,
                      # '-DPROJECT_VERSION=' + self.distribution.get_version()]

        cfg = 'Debug' if self.debug else 'Release'
        build_args = ['--config', cfg]

        if platform.system() == "Windows":
            cmake_args += ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{}={}'.format(cfg.upper(), extdir)]
            if sys.maxsize > 2**32:
                cmake_args += ['-A', 'x64']
            build_args += ['--', '/m']
        else:
            cmake_args += ['-DCMAKE_BUILD_TYPE=' + cfg]
            build_args += ['--', '-j2']

        env = os.environ.copy()
        env['CXXFLAGS'] = '{} -isystem {}'.format(env.get('CXXFLAGS', ''),
                                                  sysconfig.get_path('include'))

        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)

        subprocess.check_call(['cmake', ext.sourcedir] + cmake_args, cwd=self.build_temp, env=env)
        subprocess.check_call(['cmake', '--build', '.'] + build_args, cwd=self.build_temp)

include_dirs = os.environ.get('INCLUDE_PATH', '').split(':') + [os.path.abspath('.')]
def ecell4_cython_extension(name, libs):
    return Extension('ecell4.{}'.format(name),
                     sources = ['ecell4/{}.pyx'.format(name)],
                     include_dirs = include_dirs,
                     libraries=libs,
                     language = "c++")

cython_modules = [
    ('core', ['ecell4-core']),
    ('bd', ['ecell4-core', 'ecell4-bd']),
    ('ode', ['ecell4-core', 'ecell4-ode']),
    ('meso', ['ecell4-core', 'ecell4-meso']),
    ('egfrd', ['ecell4-core', 'ecell4-egfrd', 'greens_functions']),
    ('gillespie', ['ecell4-core', 'ecell4-gillespie']),
    ('spatiocyte', ['ecell4-core', 'ecell4-spatiocyte']),
]

shared_libs = set([])
for name, libs in cython_modules:
    for lib in libs:
        shared_libs.add(lib)
shared_libs = ['lib{}.so'.format(lib) for lib in shared_libs]

class CustomInstallCommand(install):
    def run(self):
        super(install, self).run()
        for lib in shared_libs:
            shutil.copy(os.path.join(shared_libs_dir, lib),
                        os.path.join(self.install_base, 'lib', lib))

    def get_outputs(self):
        return super(install, self).get_outputs() + \
                [os.path.join(self.install_base, 'lib', lib) for lib in shared_libs]

cython_modules = [ecell4_cython_extension(name, libs) for name, libs in cython_modules]

cwd = os.getcwd()
cython_source = os.path.abspath(os.path.join('python', 'lib'))
cython_temp = os.path.abspath(os.path.join('build', 'temp.cython'))
os.chdir(cython_source)
cython_modules = cythonize(cython_modules, build_dir=cython_temp)
os.chdir(cwd)

setup(
    name='ecell',
    version = '4.2.0',
    packages = ['ecell4'],
    package_dir = {'': cython_source},
    ext_modules=[CMakeExtension('ecell4')] + cython_modules,
    cmdclass={
        'build_ext': CMakeBuild,
        'install': CustomInstallCommand
    }
)
