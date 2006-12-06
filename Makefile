

CXXFLAGS=-Wall -O2 -g # -ftree-vectorize -ftree-vectorizer-verbose=6 -msse2
NUMPYINCLUDE=$(shell python -c "import numpy; print numpy.get_include()")

all: 
	g++ $(CXXFLAGS) -fPIC -c PlainPairGreensFunction.cpp
	g++ $(CXXFLAGS) -fPIC -I. -I/usr/include/python2.4 -I$(NUMPYINCLUDE) -c pyGFRD.cpp 
	g++ $(CXXFLAGS) -shared -fPIC PlainPairGreensFunction.o pyGFRD.o -lboost_python -lgsl -lgslcblas -lm -o gfrd.so
#	g++ -O2 -g bessel.c FullSol.c -lgsl -lgslcblas -lm -o FS.exe
#	g++ -O2 -pg -g  FullSol.c -lgsl -lgslcblas -lm -o FS.exe
#	g++ -O2 -pg -g bessel.c FullSol.c /usr/lib64/libgsl.a -lgslcblas -lm -o FS.exe


test: *.cpp
	g++ $(CXXFLAGS) -D__TEST_PLAINPAIRGREENSFUNCTION PlainPairGreensFunction.cpp -lgsl -lgslcblas -lm -o test

dist:
	( cd .. && tar -czf brown/brown.tgz brown/*.py brown/*.?pp brown/Makefile )
