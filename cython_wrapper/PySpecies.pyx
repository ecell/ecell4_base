from libcpp.vector cimport vector
from libcpp.string cimport string 

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

cdef to_PySpecies(Species *sp):
    cdef Species *new_obj = new Species("")
    r = PySpecies("")
    del r.thisptr
    r.thisptr = new_obj
    return r
