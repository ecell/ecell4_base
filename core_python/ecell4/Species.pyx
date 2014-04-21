from cython.operator cimport dereference as deref
from libcpp.string cimport string
cimport util


cdef class Species:

    def __cinit__(self, name = None):
        if name is None:
            self.thisptr = new Cpp_Species()
        else:
            self.thisptr = new Cpp_Species(<string> name)

    def __dealloc__(self):
        del self.thisptr

    def __richcmp__(Species self, Species rhs, int op):
        cdef int compare
        if deref(self.thisptr) > deref(rhs.thisptr):
            compare = 1
        elif deref(self.thisptr) < deref(rhs.thisptr):
            compare = -1
        else: # self == rhs
            compare = 0
        return util.richcmp_helper(compare, op)

    def __hash__(self):
        return hash(self.thisptr.serial())

    def serial(self):
        return self.thisptr.serial()

    def get_attribute(self, string attr_name):
        return self.thisptr.get_attribute(attr_name)

    def set_attribute(self, string name, string value):
        self.thisptr.set_attribute(name, value)

    def remove_attribute(self, string name):
        self.thisptr.remove_attribute(name)

    def has_attribute(self, string name):
        return self.thisptr.has_attribute(name)

    def add_unit(self, UnitSpecies usp):
        self.thisptr.add_unit(deref(usp.thisptr))

    def match(self, Species rhs):
        return self.thisptr.match(deref(rhs.thisptr))

    def num_units(self):
        return self.thisptr.num_units()

cdef Species Species_from_Cpp_Species(Cpp_Species *sp):
    cdef Cpp_Species *new_obj = new Cpp_Species(deref(sp))
    r = Species("")
    del r.thisptr
    r.thisptr = new_obj
    return r
