import os
import re
import sys
import glob
import platform
import sysconfig
import subprocess

from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
from setuptools.command.test import test
from distutils.version import LooseVersion


class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=''):
        super(CMakeExtension, self).__init__(name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)


class CMakeBuild(build_ext):
    def run(self):
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
            self.build_extension(ext)

    def build_extension(self, ext):
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))
        cmake_args = ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=' + extdir,
                      '-DPYTHON_EXECUTABLE=' + sys.executable]

        cfg = 'Debug' if self.debug else 'Release'
        build_args = ['--config', cfg]

        if platform.system() == "Windows":
            cmake_args += ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{}={}'.format(cfg.upper(), extdir)]
            if sys.maxsize > 2**32:
                cmake_args += ['-A', 'x64']
            build_args += ['--', '/m']
        else:
            cmake_args += ['-DCMAKE_BUILD_TYPE=' + cfg]

        if platform.system() == "Windows":
            env = os.environ.copy()
            env['CXXFLAGS'] = '{} -DVERSION_INFO=\\"{}\\" -I {}'.format(
                    env.get('CXXFLAGS', ''),
                    self.distribution.get_version(),
                    sysconfig.get_path('include'))
        else:
            env = os.environ.copy()
            env['CXXFLAGS'] = '{} -DVERSION_INFO=\\"{}\\" -isystem {}'.format(
                    env.get('CXXFLAGS', ''),
                    self.distribution.get_version(),
                    sysconfig.get_path('include'))

        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)
        subprocess.check_call(['cmake', ext.sourcedir] + cmake_args, cwd=self.build_temp, env=env)
        subprocess.check_call(['cmake', '--build', '.'] + build_args, cwd=self.build_temp)


class CustomTestCommand(test):
    def run(self):
        super(CustomTestCommand, self).run()
        build_py = self.get_finalized_command('build_ext')
        subprocess.check_call(['ctest', '--output-on-failure'], cwd=build_py.build_temp)


DESCRIPTION = (
    "A software platform for modeling, simulation and analysis of complex, "
    "heterogeneous and multi-scale systems like the cell. E-Cell has "
    "multi-algorithm, multi-timescale and multi-spatial-representation as "
    "its central feature."
)

LONG_DESCRIPTION = open("README.md").read()


setup(
    name='ecell4_base',
    version = '2.0.0b1',
    license = "the GNU General Public License v2",
    author = "Kazunari Kaizu",
    author_email = "kaizu@riken.jp",
    url = "https://github.com/ecell/ecell4-base",
    description = DESCRIPTION,
    long_description = LONG_DESCRIPTION,
    long_description_content_type='text/markdown',
    data_files = [('ecell4-licenses', glob.glob('licenses/*'))],
    ext_modules=[CMakeExtension('ecell4_base')],
    cmdclass=dict(build_ext=CMakeBuild, test=CustomTestCommand),
    zip_safe=False,
)
