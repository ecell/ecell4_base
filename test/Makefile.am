PYTHON = @PYTHON@
#PYTHONPATH = $(PYTHONPATH):$(top_srcdir)

CPP_TESTS = AllTests

PYTHON_TESTS = alltests.py

TESTS = $(CPP_TESTS) $(PYTHON_TESTS)

TESTS_ENVIRONMENT = PYTHONPATH=$(top_srcdir)


EXTRA_DIST=\
AllTests.cpp\
DynamicPriorityQueue_test.cpp\
alltests.py\
EventScheduler_test.py\
FreeGreensFunction_test.py\
FreePairGreensFunction_test.py\
FirstPassageGreensFunction_test.py\
FirstPassagePairGreensFunction_test.py

#%.py:
#	$(TESTS_ENVIRONMENT) $(PYTHON) $<

AllTests_LDFLAGS = -L$(libdir)
AllTests_LDADD   = -lboost_unit_test_framework

INCLUDES = -I$(top_srcdir)

check_PROGRAMS = $(CPP_TESTS)

AllTests_SOURCES = \
DynamicPriorityQueue_test.cpp




