PYTHON = @PYTHON@
#PYTHONPATH = $(PYTHONPATH):$(top_srcdir)

CPP_TESTS =\
DynamicPriorityQueue_test\
SphericalBesselGenerator_test

PYTHON_TESTS = alltests.py

TESTS = $(CPP_TESTS) $(PYTHON_TESTS)

TESTS_ENVIRONMENT = PYTHONPATH=$(top_srcdir)


EXTRA_DIST=\
AllTests.cpp\
DynamicPriorityQueue_test.cpp\
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
object_matrix_test.py

#%.py:
#	$(TESTS_ENVIRONMENT) $(PYTHON) $<

DynamicPriorityQueue_test_CXXFLAGS = -DBOOST_TEST_DYN_LINK
DynamicPriorityQueue_test_LDFLAGS = -L$(libdir) 
DynamicPriorityQueue_test_LDADD   = -lboost_unit_test_framework

SphericalBesselGenerator_test_CXXFLAGS = -DBOOST_TEST_DYN_LINK
SphericalBesselGenerator_test_LDFLAGS = -L$(libdir) $(GSL_LIBS)
SphericalBesselGenerator_test_LDADD   = -lboost_unit_test_framework  ../SphericalBesselGenerator.lo

INCLUDES = -I$(top_srcdir)

check_PROGRAMS = $(CPP_TESTS)

DynamicPriorityQueue_test_SOURCES = \
DynamicPriorityQueue_test.cpp

SphericalBesselGenerator_test_SOURCES = \
SphericalBesselGenerator_test.cpp




