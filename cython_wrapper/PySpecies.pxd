
from libcpp.vector cimport vector
from libcpp.string cimport string

cdef extern from "ecell4/core/Species.hpp" namespace "ecell4":
    cdef cppclass Species:
        # Constructor
        Species(string) except +
        Species(string, string)
        Species(string, string, string)
#       serial_type serial()
        string name()
        string get_attribute(string)
        void set_attribute(string,string)
        void remove_attribute(string)

cdef class PySpecies:
    cdef Species *thisptr
