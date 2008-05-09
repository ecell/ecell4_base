from numpy.distutils.core import setup, Extension

setup(
    name = 'test',
    version = '0.0.0',
    description = 'my first Cython module',
    packages = [ '' ],
    ext_modules = [
        Extension(
            'object_matrix',
            [ 'object_matrix_module.cpp', ],
            include_dirs = [ '../src', '.' ],
            #libraries = [ 'boost_python-mt', 'boost_coroutine-mt' ],
            libraries = [ 'boost_python' ],
            define_macros = [ ('USE_NUMPY', '1') ],
            language = 'c++',
            )
        ]
    )

