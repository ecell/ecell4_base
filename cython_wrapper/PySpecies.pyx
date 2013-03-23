# distutils: language = c++
# distutils: sources = ../core/Species.cpp

from libcpp.vector cimport vector
from libcpp.string cimport string 

ctypedef double Real 

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
    #cdef Species *thisptr
    def __cinit__(self, string name):
        self.thisptr = new Species(name)
    def __dealloc__(self):
        del self.thisptr 

    def name(self):
        return self.thisptr.name()
    def get_attribute(self, string attr_name):
        return self.thisptr.get_attribute(attr_name)
    def set_attribute(self, string name, string value):
        self.thisptr.set_attribute(name, value)
    def remove_attributes(self, string name):
        self.thisptr.remove_attribute(name)
