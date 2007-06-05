SUBDIRS = . test

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
AM_CXXFLAGS=@CXXFLAGS@ -Wall -g # -msse2 -mfpmath=sse  #-ffast-math 
#AM_CXXFLAGS+= -ftree-vectorize -ftree-vectorizer-verbose=5

LIBBOOSTPYTHON=-lboost_python

CFLAGS=${CXXFLAGS}

noinst_PYTHON=\
gfrdbase.py\
logger.py\
gfrdbase.py\
utils.py\
egfrd.py\
gfrd.py\
surface.py\
single.py\
dimer.py  


lib_LTLIBRARIES= _gfrd.la

noinst_HEADERS=\
Defs.hpp\
bessel.hpp\
factorial.hpp\
HalfOrderBesselGenerator.hpp\
DynamicPriorityQueue.hpp\
EventScheduler.hpp\
PyEventScheduler.hpp\
PairGreensFunction.hpp\
FirstPassageGreensFunction.hpp\
PlainPairGreensFunction.hpp\
FreePairGreensFunction.hpp\
FirstPassagePairGreensFunction.hpp



_gfrd_la_SOURCES=\
PlainPairGreensFunction.cpp\
FreePairGreensFunction.cpp\
FirstPassageGreensFunction.cpp\
FirstPassagePairGreensFunction.cpp\
bessel.cpp\
pyGFRD.cpp



_gfrd_la_LDFLAGS=-module -export-dynamic -avoid-version $(no_undefined)

#GSL_LIBS=/usr/lib/libgsl.a /usr/lib/libgslcblas.a

_gfrd_la_LIBADD=$(LIBBOOSTPYTHON) $(LIBPYTHON) $(GSL_LIBS)

GFRD_SO=_gfrd$(LTDL_SHLIB_EXT)

all-local: $(GFRD_SO)

$(GFRD_SO): _gfrd.la 
	ln -sf .libs/$(GFRD_SO) ./$(GFRD_SO)

clean-local:
	rm -f $(GFRD_SO)


PYCHECKER = @PYCHECKER@
PYCHECKER_FLAGS = --limit 0 

pycheck:
	$(PYCHECKER) $(PYCHECKER_FLAGS) $(top_srcdir)/egfrd.py


EXTRA_DIST=\
autogen.sh

