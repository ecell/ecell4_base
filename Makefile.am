ACLOCAL_AMFLAGS=-I m4

SUBDIRS = doc test

TOPSRCDIR = ${top_srcdir}

LTDL_SHLIB_EXT = @LTDL_SHLIB_EXT@

PYTHON = @PYTHON@
PYTHON_PREFIX = @PYTHON_PREFIX@
PYTHON_VERSION = @PYTHON_VERSION@

NUMPY_INCLUDE_DIR=@NUMPY_INCLUDE_DIR@
PYTHON_INCLUDES = @PYTHON_INCLUDES@

pkgpyexecdir = @pkgpyexecdir@

LIBPYTHON = -lpython$(PYTHON_VERSION)

INCLUDES = ${PYTHON_INCLUDES} -I${NUMPY_INCLUDE_DIR}
AM_CXXFLAGS = @CXXFLAGS@ @GSL_CFLAGS@ -Wall -g
LIBBOOSTPYTHON = -lboost_python

CFLAGS = @GSL_CFLAGS@

PYCHECKER = @PYCHECKER@
PYCHECKER_FLAGS = --limit 0 

noinst_PYTHON = \
	bd.py\
	cObjectMatrix.py\
	gfrdbase.py\
	logger.py\
	gfrdbase.py\
	utils.py\
	egfrd.py\
	gfrd.py\
	surface.py\
	single.py\
	dimer.py  


lib_LTLIBRARIES = _gfrd.la object_matrix.la

noinst_HEADERS = \
	BasicPairGreensFunction.hpp\
	Defs.hpp\
	DynamicPriorityQueue.hpp\
	EventScheduler.hpp\
	FirstPassageGreensFunction.hpp\
	FirstPassageNoCollisionPairGreensFunction.hpp\
	FirstPassagePairGreensFunction.hpp\
	FreeGreensFunction.hpp\
	FreePairGreensFunction.hpp\
	GreensFunction.hpp\
	HalfOrderBesselGenerator.hpp\
	PairGreensFunction.hpp\
	PyEventScheduler.hpp\
	SphericalBesselGenerator.hpp\
	array_cast.hpp\
	array_helper.hpp\
	bessel.hpp\
	distance.hpp\
	factorial.hpp\
	filters.hpp\
	findRoot.hpp\
	freeFunctions.hpp\
	funcSum.hpp\
	get_mapper_mf.hpp\
	object_container.hpp\
	position.hpp\
	sphere.hpp\
	utils.hpp\
	vmap.hpp\
	wrapped_multi_array.hpp\
	peer/ObjectContainer.hpp\
	peer/enumerator.hpp\
	peer/generator_support.hpp\
	peer/tuple_converters.hpp\
	peer/utils.hpp\
	peer/numpy/ndarray_converters.hpp\
	peer/numpy/pyarray_backed_allocator.hpp\
	peer/numpy/type_mappings.hpp

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


_gfrd_la_LDFLAGS = -module -export-dynamic -avoid-version $(no_undefined)
_gfrd_la_LIBADD = $(LIBBOOSTPYTHON) $(LIBPYTHON) $(GSL_LIBS)

object_matrix_la_LDFLAGS = -module -export-dynamic -avoid-version $(no_undefined)
object_matrix_la_LIBADD = $(LIBBOOSTPYTHON) $(LIBPYTHON) $(GSL_LIBS)
object_matrix_la_SOURCES = object_matrix_module.cpp


nodist__gfrd_la_SOURCES = \
	SphericalBesselTable.hpp

BUILT_SOURCES = \
	SphericalBesselTable.hpp

CLEANFILES = \
	SphericalBesselTable.hpp

SphericalBesselTable.hpp: make_sjy_table.py
	$(PYTHON) make_sjy_table.py SphericalBesselTable.hpp


object_matrix.so: object_matrix.la
	ln -sf .libs/object_matrix.so ./object_matrix.so

_gfrd.so: _gfrd.la 
	ln -sf .libs/_gfrd.so ./_gfrd.so

all-local: _gfrd.so object_matrix.so

clean-local:
	rm -f _gfrd.so


pycheck:
	$(PYCHECKER) $(PYCHECKER_FLAGS) $(top_srcdir)/egfrd.py


dist-hook:
	rm -rf `find $(distdir)/samples -name '.svn'`


EXTRA_DIST = \
	autogen.sh\
	samples\
	make_sjy_table.py
