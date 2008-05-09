
LTDL_SHLIB_EXT=@LTDL_SHLIB_EXT@

PYTHON=@PYTHON@
#PYTHON_PREFIX=@PYTHON_PREFIX@
PYTHON_VERSION=@PYTHON_VERSION@

NUMPY_INCLUDE_DIR=@NUMPY_INCLUDE_DIR@
PYTHON_INCLUDES=@PYTHON_INCLUDES@

pkgpyexecdir=@pkgpyexecdir@

LIBPYTHON=-lpython$(PYTHON_VERSION)

INCLUDES=${PYTHON_INCLUDES} -I${NUMPY_INCLUDE_DIR} -I../src
AM_CXXFLAGS=@CXXFLAGS@ -Wall -g # -msse2 -mfpmath=sse  #-ffast-math 


LIBBOOSTPYTHON=-lboost_python

lib_LTLIBRARIES= object_matrix.la

object_matrix_la_SOURCES=\
object_matrix_module.cpp

object_matrix_la_LDFLAGS=-module -export-dynamic -avoid-version $(no_undefined)
object_matrix_la_LIBADD=$(LIBBOOSTPYTHON) $(LIBPYTHON)
OBJECT_MATRIX_SO=object_matrix$(LTDL_SHLIB_EXT)



$(OBJECT_MATRIX_SO): object_matrix.la 
	ln -sf .libs/$(OBJECT_MATRIX_SO) ./$(OBJECT_MATRIX_SO)

all-local: $(OBJECT_MATRIX_SO)

clean-local:
	rm -f  $(OBJECT_MATRIX_SO)

