PYTHON = @PYTHON@
#PYTHONPATH = $(PYTHONPATH):$(top_srcdir)

CPP_TESTS =\
DynamicPriorityQueue_test\
SphericalBesselGenerator_test\
array_helper_test\
filters_test\
object_container_test\
vmap_test

PYTHON_TESTS = alltests.py

TESTS = $(CPP_TESTS) $(PYTHON_TESTS)

TESTS_ENVIRONMENT = PYTHONPATH=$(top_srcdir)


EXTRA_DIST=\
AllTests.cpp\
DynamicPriorityQueue_test.cpp\
array_helper_test.cpp\
filters_test.cpp\
object_container_test.cpp\
vmap_test.cpp\
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

object_container_test_SOURCES = object_container_test.cpp

vmap_test_SOURCES = vmap_test.cpp

