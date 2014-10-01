noinst_PROGRAMS = hardbody

INCLUDES = -I$(top_srcdir) $(PYTHON_INCLUDES)
AM_CXXFLAGS = @CXXFLAGS@ @GSL_CFLAGS@

hardbody_SOURCES = hardbody.cpp ../../NetworkRules.cpp ../../BasicNetworkRulesImpl.cpp ../../Logger.cpp ../../freeFunctions.cpp

hardbody_LDADD = $(GSL_LIBS)
