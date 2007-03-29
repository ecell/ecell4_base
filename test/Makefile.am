PYCHECKER = @PYCHECKER@


TESTS = AllTests 

#pycheck


INCLUDES = -I$(top_srcdir)

check_PROGRAMS = AllTests

AllTests_SOURCES = AllTests.cpp\
DynamicPriorityQueue_test.cpp


AllTests_LDFLAGS = -L$(libdir)
AllTests_LDADD   = $(CPPUNIT_LIBS)


PYCHECKER_FLAGS = --limit 0

pycheck:
	$(PYCHECKER) $(PYCHECKER_FLAGS) $(top_srcdir)/egfrd.py
