from cython.operator cimport dereference as deref
from libcpp.string cimport string
cimport util


cdef class UnitSpecies:

    def __cinit__(self, string name):
        self.thisptr = new Cpp_UnitSpecies(name)

    def __dealloc__(self):
        del self.thisptr

    def __richcmp__(UnitSpecies self, UnitSpecies rhs, int op):
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

    def name(self):
        return self.thisptr.name()

    def serial(self):
        return self.thisptr.serial()

    def add_site(self, string name, string state, string bond):
        return self.thisptr.add_site(name, state, bond)

cdef UnitSpecies UnitSpecies_from_Cpp_UnitSpecies(Cpp_UnitSpecies *sp):
    cdef Cpp_UnitSpecies *new_obj = new Cpp_UnitSpecies(deref(sp))
    r = UnitSpecies("")
    del r.thisptr
    r.thisptr = new_obj
    return r
