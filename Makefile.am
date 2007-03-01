TOPSRCDIR=${top_srcdir}

LTDL_SHLIB_EXT=@LTDL_SHLIB_EXT@

PYTHON=@PYTHON@
#PYTHON_PREFIX=@PYTHON_PREFIX@
PYTHON_VERSION=@PYTHON_VERSION@

NUMPY_INCLUDE_DIR=$(shell python -c "import numpy; print numpy.get_include()")
PYTHON_INCLUDES=@PYTHON_INCLUDES@

pkgpyexecdir=@pkgpyexecdir@

LIBPYTHON=-lpython$(PYTHON_VERSION)

INCLUDES=${PYTHON_INCLUDES} -I${NUMPY_INCLUDE_DIR}
AM_CXXFLAGS=@CXXFLAGS@ -g # -ffast-math # -ftree-vectorize -ftree-vectorizer-verbose=5

LIBBOOSTPYTHON=-lboost_python

CFLAGS=${CXXFLAGS}

#pkgpyexec_LTLIBRARIES = _ecs.la _emc.la


noinst_LTLIBRARIES= bessel.la 
lib_LTLIBRARIES= _gfrd.la

noinst_HEADERS=\
Defs.hpp\
bessel.hpp\
HalfOrderBesselGenerator.hpp\
DynamicPriorityQueue.hpp\
EventScheduler.hpp\
PairGreensFunction.hpp\
FirstPassageGreensFunction.hpp\
PlainPairGreensFunction.hpp\
FirstPassagePairGreensFunction.hpp


bessel_la_SOURCES=\
bessel.cpp
bessel_la_LDFLAGS=-module -export-dynamic -avoid-version $(no_undefined)


_gfrd_la_SOURCES=\
PlainPairGreensFunction.cpp\
FirstPassageGreensFunction.cpp\
FirstPassagePairGreensFunction.cpp\
pyGFRD.cpp

_gfrd_la_LDFLAGS=-module -export-dynamic -avoid-version $(no_undefined)

#GSL_LIBS=/usr/lib64/libgsl.a /usr/lib64/libgslcblas.a

_gfrd_la_LIBADD=$(LIBBOOSTPYTHON) $(LIBPYTHON) $(GSL_LIBS) bessel.la

_gfrd_la_DEPENDENCIES=bessel.la

GFRD_SO=_gfrd$(LTDL_SHLIB_EXT)

all-local: $(GFRD_SO)

$(GFRD_SO): _gfrd.la 
	ln -sf .libs/$(GFRD_SO) ./$(GFRD_SO)

clean-local:
	rm -f $(GFRD_SO)

