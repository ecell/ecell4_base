ACLOCAL_AMFLAGS=-I m4

SUBDIRS = binding doc test 

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
AM_CXXFLAGS = @CXXFLAGS@ @GSL_CFLAGS@
LIBBOOSTPYTHON = -l@BOOST_PYTHON_LIBNAME@
LIBBOOSTREGEX = -l@BOOST_REGEX_LIBNAME@

CFLAGS = @GSL_CFLAGS@ @BOOST_CPPFLAGS@

PYCHECKER = @PYCHECKER@
PYCHECKER_FLAGS = --limit 0 

noinst_PYTHON = \
	bd.py\
	gfrdbase.py\
	logger.py\
	gfrdbase.py\
	utils.py\
	egfrd.py\
	gfrd.py\
	surface.py\
	single.py\
	dimer.py  


lib_LTLIBRARIES = _gfrd.la _greens_functions.la

noinst_HEADERS = \
	abstract_set.hpp\
	BasicNetworkRulesImpl.hpp\
	GreensFunction3DRadInf.hpp\
	BDPropagator.hpp\
	BDSimulator.hpp\
	bessel.hpp\
	Box.hpp\
	ConsoleAppender.hpp\
	Cylinder.hpp\
	Defs.hpp\
	DomainFactory.hpp\
	Domain.hpp\
	DomainID.hpp\
	DynamicPriorityQueue.hpp\
	EGFRDSimulator.hpp\
	EventScheduler.hpp\
	exceptions.hpp\
	factorial.hpp\
	filters.hpp\
	findRoot.hpp\
	GreensFunction1DRadAbs.hpp\
	GreensFunction1DAbsAbs.hpp\
	GreensFunction3DAbsSym.hpp\
	GreensFunction3DAbs.hpp\
	GreensFunction3DRadAbs.hpp\
	GreensFunction3DRadAbsBase.hpp\
	freeFunctions.hpp\
	GreensFunction3DSym.hpp\
	GreensFunction3D.hpp\
	funcSum.hpp\
	generator.hpp\
	geometry.hpp\
	GreensFunction.hpp\
	GSLRandomNumberGenerator.hpp\
	gsl_rng_base.hpp\
	HalfOrderBesselGenerator.hpp\
	Identifier.hpp\
	linear_algebra.hpp\
	Logger.hpp\
	MatrixSpace.hpp\
	Model.hpp\
	Multi.hpp\
	NetworkRules.hpp\
	NetworkRulesWrapper.hpp\
	PairGreensFunction.hpp\
	Pair.hpp\
	ParticleContainerBase.hpp\
	ParticleContainer.hpp\
	Particle.hpp\
	ParticleID.hpp\
	ParticleSimulator.hpp\
	Point.hpp\
	PyEventScheduler.hpp\
	ReactionRule.hpp\
	ReactionRuleInfo.hpp\
	Region.hpp\
	SerialIDGenerator.hpp\
	Shape.hpp\
	Shell.hpp\
	ShellID.hpp\
	Single.hpp\
	sorted_list.hpp\
	SpeciesInfo.hpp\
	SpeciesType.cpp\
	SpeciesType.hpp\
	SpeciesTypeID.hpp\
	Sphere.hpp\
	SphericalBesselGenerator.hpp\
	SphericalBesselTable.hpp\
	Structure.hpp\
	StructureUtils.hpp\
	StructureType.hpp\
	Surface.hpp\
	Transaction.hpp\
	twofold_container.hpp\
	Vector3.hpp\
	World.hpp\
	peer/numpy/scalar_converters.hpp \
	peer/numpy/ndarray_converters.hpp \
	peer/numpy/pyarray_backed_allocator.hpp \
	peer/numpy/type_mappings.hpp \
	peer/numpy/wrapped_multi_array.hpp \
	peer/set_indexing_suite.hpp \
	peer/compat.h \
	peer/pickle_support.hpp \
	peer/converters/sequence.hpp \
	peer/converters/tuple.hpp \
	peer/converters/tuple/to_python.hpp \
	peer/converters/tuple/from_python.hpp \
	peer/converters/iterator.hpp \
	peer/converters/generator/to_python.hpp \
	peer/converters/generator/from_python.hpp \
	peer/converters/iterator/to_python.hpp \
	peer/converters/sequence/to_python.hpp \
	peer/converters/sequence/from_python.hpp \
	peer/wrappers/range/pyiterable_range.hpp \
	peer/wrappers/range/stl_container_wrapper.hpp \
	peer/wrappers/exception/exception_wrapper.hpp \
	peer/wrappers/generator/pyiterator_generator.hpp \
	peer/wrappers/generator/generator_wrapper.hpp \
	peer/wrappers/iterator/pyseq_iterator.hpp \
	peer/wrappers/iterator/stl_iterator_wrapper.hpp \
	peer/py_hash_support.hpp \
	peer/utils.hpp \
	utils/array_helper.hpp\
	utils/array_traits.hpp\
	utils/fun_composition.hpp\
	utils/fun_wrappers.hpp\
	utils/get_default_impl.hpp\
	utils/get_mapper_mf.hpp\
	utils.hpp\
	utils/memberwise_compare.hpp\
	utils/pair.hpp\
	utils/pointer_preds.hpp\
	utils/range.hpp\
	utils/range_support.hpp\
	utils/reference_or_instance.hpp\
	utils/unassignable_adapter.hpp

_gfrd_la_CPPFLAGS = -DPY_ARRAY_UNIQUE_SYMBOL=PyArray_API

_gfrd_la_SOURCES=\
	BasicNetworkRulesImpl.cpp\
	findRoot.cpp\
	freeFunctions.cpp\
	funcSum.cpp\
	GreensFunction3DAbs.cpp\
	GreensFunction3DAbsSym.cpp\
	GreensFunction3D.cpp\
	GreensFunction3DRadAbsBase.cpp\
	GreensFunction3DRadAbs.cpp\
	GreensFunction3DRadInf.cpp\
	GreensFunction3DSym.cpp\
	Logger.cpp\
	ConsoleAppender.cpp\
	Model.cpp\
	NetworkRules.cpp\
	ParticleModel.cpp\
	pyGFRD.cpp\
	SpeciesType.cpp\
	SphericalBesselGenerator.cpp\
	StructureType.cpp\
	utils.cpp

_greens_functions_la_SOURCES=\
	findRoot.cpp\
	funcSum.cpp\
	greens_functions.cpp\
	SphericalBesselGenerator.cpp\
	freeFunctions.cpp\
	GreensFunction1DAbsAbs.cpp\
	GreensFunction1DRadAbs.cpp\
	GreensFunction3DSym.cpp\
	GreensFunction3DAbsSym.cpp\
	GreensFunction3DRadInf.cpp\
	GreensFunction3D.cpp\
	GreensFunction3DRadAbsBase.cpp\
	GreensFunction3DRadAbs.cpp\
	GreensFunction3DAbs.cpp\
	Logger.cpp\
	ConsoleAppender.cpp

_gfrd_la_LDFLAGS = -module -export-dynamic -avoid-version $(no_undefined)
_gfrd_la_LIBADD = binding/libbinding_utils.la
_gfrd_la_LIBADD += $(LIBBOOSTPYTHON) $(LIBPYTHON) $(GSL_LIBS)

_greens_functions_la_LDFLAGS = -module -export-dynamic -avoid-version $(no_undefined)
_greens_functions_la_LIBADD = $(LIBBOOSTPYTHON) $(LIBPYTHON) $(GSL_LIBS)

#nodist__gfrd_la_SOURCES = \
#	SphericalBesselTable.hpp

nodist__greens_functions_la_SOURCES = \
	SphericalBesselTable.hpp\
	CylindricalBesselTable.hpp

BUILT_SOURCES = \
	SphericalBesselTable.hpp\
	CylindricalBesselTable.hpp

CLEANFILES = \
	SphericalBesselTable.hpp\
	CylindricalBesselTable.hpp

SphericalBesselTable.hpp: make_sjy_table.py
	$(PYTHON) make_sjy_table.py SphericalBesselTable.hpp

CylindricalBesselTable.hpp: make_cjy_table.py
	$(PYTHON) make_cjy_table.py CylindricalBesselTable.hpp


_gfrd.so: _gfrd.la 
	ln -sf .libs/_gfrd.so ./_gfrd.so

_greens_functions.so: _greens_functions.la
	ln -sf .libs/_greens_functions.so ./_greens_functions.so

all-local: _gfrd.so _greens_functions.so

clean-local:
	rm -f _gfrd.so _greens_functions.so


pycheck:
	$(PYCHECKER) $(PYCHECKER_FLAGS) $(top_srcdir)/egfrd.py


dist-hook:
	rm -rf `find $(distdir)/samples -name '.svn'`


EXTRA_DIST = \
	autogen.sh\
	samples\
	make_sjy_table.py
