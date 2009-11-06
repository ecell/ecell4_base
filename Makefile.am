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
LIBBOOSTPYTHON = -l@BOOST_PYTHON_LIBNAME@

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


lib_LTLIBRARIES = _gfrd.la

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
	array_traits.hpp\
	array_helper.hpp\
	bessel.hpp\
	factorial.hpp\
	filters.hpp\
	findRoot.hpp\
	freeFunctions.hpp\
	funcSum.hpp\
	get_mapper_mf.hpp\
	MatrixSpace.hpp\
	Vector3.hpp\
	Sphere.hpp\
	NetworkRules.hpp\
	ReactionRule.hpp\
	BasicNetworkRulesImpl.hpp\
	Identifier.hpp\
	Model.hpp\
	utils.hpp\
	vmap.hpp\
	peer/Exception.hpp\
	peer/MatrixSpace.hpp\
	peer/GeneratorIteratorWrapper.hpp\
	peer/tuple_converters.hpp\
	peer/utils.hpp\
	peer/numpy/wrapped_multi_array.hpp\
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
	pyGFRD.cpp\
	Model.cpp\
	BasicNetworkRulesImpl.cpp\
	NetworkRules.cpp\
	utils.cpp


_gfrd_la_LDFLAGS = -module -export-dynamic -avoid-version $(no_undefined)
_gfrd_la_LIBADD = $(LIBBOOSTPYTHON) $(LIBPYTHON) $(GSL_LIBS)

nodist__gfrd_la_SOURCES = \
	SphericalBesselTable.hpp

BUILT_SOURCES = \
	SphericalBesselTable.hpp

CLEANFILES = \
	SphericalBesselTable.hpp

SphericalBesselTable.hpp: make_sjy_table.py
	$(PYTHON) make_sjy_table.py SphericalBesselTable.hpp


_gfrd.so: _gfrd.la 
	ln -sf .libs/_gfrd.so ./_gfrd.so

all-local: _gfrd.so

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
