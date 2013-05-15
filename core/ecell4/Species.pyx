from cython.operator cimport dereference as deref
from libcpp.string cimport string


cdef class Species:

    def __cinit__(self, string name):
        self.thisptr = new Cpp_Species(name)

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

cdef Species Cpp_Species_to_Species(Cpp_Species *sp):
    cdef Cpp_Species *new_obj = new Cpp_Species(deref(sp))
    r = Species("")
    del r.thisptr
    r.thisptr = new_obj
    return r
