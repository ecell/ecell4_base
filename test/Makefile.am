PYTHON = @PYTHON@
#PYTHONPATH = $(PYTHONPATH):$(top_srcdir)

CPP_TESTS =\
DynamicPriorityQueue_test\
SphericalBesselGenerator_test\
array_helper_test\
filters_test\
MatrixSpace_test\
MatrixSpaceWithCylinders_test\
World_test\
model_test\
Vector3_test\
BDPropagator_test\
range_support_test\
py_range_converters_test \
linear_algebra_test\
geometry_test\
StructureUtils_test\
sorted_list_test\
pointer_as_ref_test\
EGFRDSimulator_test

PYTHON_TESTS = \
	BDSimulator_test.py \
	CylindricalShellContainer_test.py \
	EGFRDSimulator_test.py \
	EventScheduler_test.py \
	GreensFunction3DRadAbs_test.py \
	GreensFunction3DRadInf_test.py \
	GreensFunction3DAbs_test.py \
	GreensFunction3DAbsSym_test.py \
	GreensFunction3DSym_test.py \
	GreensFunction3D_test.py \
	Model_test.py \
	NetworkRules_test.py \
	NetworkRulesWrapper_test.py \
	ReactionRule_test.py \
	SphericalShellContainer_test.py \
	freeFunctions_test.py \
	CylindricalSurface_test.py \
	PlanarSurface_test.py \
	utils_test.py \
	ReactionRecord_test.py

#GreensFunction1DAbsAbs_test.py \
#GreensFunction1DRadAbs_test.py

TESTS = $(CPP_TESTS) $(PYTHON_TESTS)

TESTS_ENVIRONMENT = PYTHONPATH=$(top_srcdir)


EXTRA_DIST=\
AllTests.cpp\
DynamicPriorityQueue_test.cpp\
array_helper_test.cpp\
filters_test.cpp\
MatrixSpace_test.cpp\
alltests.py\
utils_test.py\
freeFunctions_test.py\
EventScheduler_test.py\
GreensFunction1DAbsAbs_test.py\
GreensFunction1DRadAbs_test.py\
GreensFunction3DSym_test.py\
GreensFunction3D_test.py\
GreensFunction3DRadInf_test.py\
GreensFunction3DAbsSym_test.py\
GreensFunction3DRadAbs_test.py\
EGFRDSimulator_test.py\
SphericalShellContainer_test.py\
CylindricalShellContainer_test.py\
CylindricalSurface_test.py \
PlanarSurface_test.py \
Model_test.py\
NetworkRules_test.py\
ReactionRule_test.py\
ReactionRecord_test.py

#%.py:
#	$(TESTS_ENVIRONMENT) $(PYTHON) $<

INCLUDES = -I$(top_srcdir) $(PYTHON_INCLUDES)
AM_CXXFLAGS = @CXXFLAGS@ @GSL_CFLAGS@ -Wall -g

check_PROGRAMS = $(CPP_TESTS)

DynamicPriorityQueue_test_SOURCES = \
DynamicPriorityQueue_test.cpp

SphericalBesselGenerator_test_LDADD = $(GSL_LIBS)
SphericalBesselGenerator_test_SOURCES = \
SphericalBesselGenerator_test.cpp ../SphericalBesselGenerator.cpp

array_helper_test_SOURCES = array_helper_test.cpp

filters_test_SOURCES = filters_test.cpp

MatrixSpace_test_SOURCES = MatrixSpace_test.cpp
MatrixSpace_test_LDADD = $(GSL_LIBS)

MatrixSpaceWithCylinders_test_SOURCES = MatrixSpaceWithCylinders_test.cpp

World_test_SOURCES = World_test.cpp

model_test_SOURCES = model_test.cpp ../Model.cpp ../NetworkRules.cpp ../BasicNetworkRulesImpl.cpp ../SpeciesType.cpp

Vector3_test_SOURCES = Vector3_test.cpp ../Vector3.hpp

BDPropagator_test_SOURCES = BDPropagator_test.cpp ../BasicNetworkRulesImpl.cpp ../NetworkRules.cpp ../Logger.cpp ../ConsoleAppender.cpp ../freeFunctions.cpp ../BDPropagator.hpp
BDPropagator_test_LDADD = $(GSL_LIBS)

range_support_test_SOURCES = range_support_test.cpp ../utils/range.hpp ../utils/range_support.hpp

linear_algebra_test_SOURCES = linear_algebra_test.cpp ../linear_algebra.hpp

geometry_test_SOURCES = geometry_test.cpp ../geometry.hpp

StructureUtils_test_SOURCES = StructureUtils_test.cpp ../StructureUtils.hpp
StructureUtils_test_LDADD = $(GSL_LIBS)
py_range_converters_test_SOURCES = py_range_converters_test.cpp ../peer/range_converters.hpp 

pointer_as_ref_test_SOURCES = pointer_as_ref_test.cpp ../utils/pointer_as_ref.hpp

EGFRDSimulator_test_SOURCES = EGFRDSimulator_test.cpp ../EGFRDSimulator.hpp ../Model.cpp ../NetworkRules.cpp ../BasicNetworkRulesImpl.cpp ../SpeciesType.cpp ../freeFunctions.cpp ../Logger.cpp ../ConsoleAppender.cpp ../GreensFunction3D.cpp ../GreensFunction3DAbs.cpp ../GreensFunction3DAbsSym.cpp ../GreensFunction3DRadAbs.cpp ../GreensFunction3DRadAbsBase.cpp ../GreensFunction3DRadInf.cpp ../GreensFunction3DSym.cpp ../SphericalBesselGenerator.cpp ../CylindricalBesselGenerator.cpp ../funcSum.cpp ../findRoot.cpp ../ParticleModel.cpp ../StructureType.cpp
EGFRDSimulator_test_LIBS = -l@BOOST_REGEX_LIBNAME@ -l@BOOST_DATE_TIME_LIBNAME@
EGFRDSimulator_test_CPPFLAGS = -DDEBUG

EGFRDSimulator_test_LDADD = $(GSL_LIBS)

sorted_list_test_SOURCES = sorted_list_test.cpp ../sorted_list.hpp
sorted_list_test_LDADD = $(GSL_LIBS)

