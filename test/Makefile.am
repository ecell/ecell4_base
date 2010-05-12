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
range_support_test

PYTHON_TESTS = \
	BDSimulator_test.py \
	BasicPairGreensFunction_test.py \
	CylindricalShellContainer_test.py \
	EGFRDSimulator_test.py \
	EventScheduler_test.py \
	FirstPassageNoCollisionPairGreensFunction_test.py \
	FirstPassagePairGreensFunction_test.py \
	FirstPassageGreensFunction_test.py \
	FreeGreensFunction_test.py \
	FreePairGreensFunction_test.py \
	Model_test.py \
	NetworkRules_test.py \
	NetworkRulesWrapper_test.py \
	ReactionRule_test.py \
	SphericalShellContainer_test.py \
	freeFunctions_test.py \
	utils_test.py

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
FreeGreensFunction_test.py\
FreePairGreensFunction_test.py\
BasicPairGreensFunction_test.py\
FirstPassageGreensFunction_test.py\
FirstPassagePairGreensFunction_test.py\
EGFRDSimulator_test.py\
SphericalShellContainer_test.py\
CylindricalShellContainer_test.py\
Model_test.py\
NetworkRules_test.py\
ReactionRule_test.py

#%.py:
#	$(TESTS_ENVIRONMENT) $(PYTHON) $<

INCLUDES = -I$(top_srcdir)
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

MatrixSpaceWithCylinders_test_SOURCES = MatrixSpaceWithCylinders_test.cpp

World_test_SOURCES = World_test.cpp

model_test_SOURCES = model_test.cpp ../Model.cpp ../NetworkRules.cpp ../BasicNetworkRulesImpl.cpp

Vector3_test_SOURCES = Vector3_test.cpp ../Vector3.hpp

BDPropagator_test_SOURCES = BDPropagator_test.cpp ../BasicNetworkRulesImpl.cpp ../NetworkRules.cpp ../Logger.cpp ../freeFunctions.cpp ../BDPropagator.hpp
BDPropagator_test_LDADD = $(GSL_LIBS)

range_support_test_SOURCES = range_support_test.cpp ../utils/range.hpp ../utils/range_support.hpp

