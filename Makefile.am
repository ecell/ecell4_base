SUBDIRS = doc test object_matrix

TOPSRCDIR=${top_srcdir}

LTDL_SHLIB_EXT=@LTDL_SHLIB_EXT@

PYTHON=@PYTHON@
#PYTHON_PREFIX=@PYTHON_PREFIX@
PYTHON_VERSION=@PYTHON_VERSION@

NUMPY_INCLUDE_DIR=@NUMPY_INCLUDE_DIR@
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
funcSum.hpp\
findRoot.hpp\	
freeFunctions.hpp\
factorial.hpp\
SphericalBesselGenerator.hpp\
DynamicPriorityQueue.hpp\
EventScheduler.hpp\
PyEventScheduler.hpp\
GreensFunction.hpp\
PairGreensFunction.hpp\
FreeGreensFunction.hpp\
FirstPassageGreensFunction.hpp\
BasicPairGreensFunction.hpp\
FreePairGreensFunction.hpp\
FirstPassagePairGreensFunction.hpp\
FirstPassageNoCollisionPairGreensFunction.hpp



_gfrd_la_SOURCES=\
funcSum.cpp\
findRoot.cpp\
freeFunctions.cpp\
SphericalBesselGenerator.cpp\
FreeGreensFunction.cpp\
FirstPassageGreensFunction.cpp\
BasicPairGreensFunction.cpp\
FreePairGreensFunction.cpp\
FirstPassagePairGreensFunction.cpp\
FirstPassageNoCollisionPairGreensFunction.cpp\
pyGFRD.cpp

_gfrd_la_LDFLAGS=-module -export-dynamic -avoid-version $(no_undefined)
_gfrd_la_LIBADD=$(LIBBOOSTPYTHON) $(LIBPYTHON) $(GSL_LIBS)
GFRD_SO=_gfrd$(LTDL_SHLIB_EXT)


OBJECT_MATRIX_SO=object_matrix.so
$(OBJECT_MATRIX_SO):
	ln -sf object_matrix/boost.python/$(OBJECT_MATRIX_SO) ./$(OBJECT_MATRIX_SO)


all-local: $(GFRD_SO) $(OBJECT_MATRIX_SO)

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

