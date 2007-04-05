PYTHON = @PYTHON@
#PYTHONPATH = $(PYTHONPATH):$(top_srcdir)

CPP_TESTS = AllTests

PYTHON_TESTS = alltests.py

TESTS = $(CPP_TESTS) $(PYTHON_TESTS)

TESTS_ENVIRONMENT = PYTHONPATH=$(top_srcdir)

#%.py:
#	$(TESTS_ENVIRONMENT) $(PYTHON) $<


INCLUDES = -I$(top_srcdir)

check_PROGRAMS = $(CPP_TESTS)

AllTests_SOURCES = AllTests.cpp\
DynamicPriorityQueue_test.cpp


AllTests_LDFLAGS = -L$(libdir)
AllTests_LDADD   = $(CPPUNIT_LIBS)

