PYTHON = @PYTHON@
#PYTHONPATH = $(PYTHONPATH):$(top_srcdir)

CPP_TESTS =\
DynamicPriorityQueue_test\
SphericalBesselTable_test

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
_object_matrix_test.py

#%.py:
#	$(TESTS_ENVIRONMENT) $(PYTHON) $<

DynamicPriorityQueue_test_CXXFLAGS = -DBOOST_TEST_DYN_LINK
DynamicPriorityQueue_test_LDFLAGS = -L$(libdir) 
DynamicPriorityQueue_test_LDADD   = -lboost_unit_test_framework

SphericalBesselTable_test_CXXFLAGS = -DBOOST_TEST_DYN_LINK
SphericalBesselTable_test_LDFLAGS = -L$(libdir) $(GSL_LIBS)
SphericalBesselTable_test_LDADD   = -lboost_unit_test_framework

INCLUDES = -I$(top_srcdir)

check_PROGRAMS = $(CPP_TESTS)

DynamicPriorityQueue_test_SOURCES = \
DynamicPriorityQueue_test.cpp

SphericalBesselTable_test_SOURCES = \
SphericalBesselTable_test.cpp




