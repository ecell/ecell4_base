SRCDIR=../src
CPPFLAGS=-I$(SRCDIR)
CXXFLAGS=-g -Werror
LDFLAGS=

all: vmap_test object_container_test array_helper_test filters_test position_test

clean:
	rm -f *.o *.so *.lo *.loT
	rm -f vmap_test object_container_test array_helper_test filters_test \
	      position_test

$(SRCDIR)/object_container.hpp: $(SRCDIR)/vmap.hpp $(SRCDIR)/sphere.hpp $(SRCDIR)/array_helper.hpp

$(SRCDIR)/sphere.hpp: $(SRCDIR)/position.hpp

vmap_test: vmap_test.cpp $(SRCDIR)/vmap.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ vmap_test.cpp

object_container_test: object_container_test.cpp $(SRCDIR)/object_container.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ object_container_test.cpp

array_helper_test: array_helper_test.cpp $(SRCDIR)/array_helper.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ array_helper_test.cpp

filters_test: filters_test.cpp $(SRCDIR)/filters.hpp $(SRCDIR)/object_container.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ filters_test.cpp

position_test: position_test.cpp $(SRCDIR)/position.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ position_test.cpp

.PHONY: all clean
