from distutils.core import setup, Extension
#from Pyrex.Distutils.build_ext import build_ext
from Cython.Distutils.build_ext import build_ext

setup(
    name = 'test',
    version = '0.0.0',
    description = 'my first Cython module',
    packages = [ '' ],
    ext_modules = [
        Extension(
            'object_matrix',
            [ '__init__.pyx', ],
            include_dirs = [ '../src' ],
            libraries = [ 'boost_coroutine-mt' ],
            define_macros = [ ('USE_NUMPY', '1') ],
            language = 'c++',
            )
        ],
    cmdclass = { 'build_ext': build_ext },
    )

