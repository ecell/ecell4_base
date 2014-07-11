from cython.operator cimport dereference as deref, preincrement as inc
from libcpp.string cimport string
from cython cimport address
cimport util

cimport context


cdef class Species:

    def __cinit__(self, serial=None, radius=None, D=None):
        if serial is None:
            self.thisptr = new Cpp_Species()
        elif radius is None or D is None:
            self.thisptr = new Cpp_Species(<string> serial)
        else:
            self.thisptr = new Cpp_Species(
                <string> serial, <string> radius, <string> D)

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

    def list_attributes(self):
        return self.thisptr.list_attributes()

    def add_unit(self, UnitSpecies usp):
        self.thisptr.add_unit(deref(usp.thisptr))

    def match(self, Species rhs):
        return self.thisptr.match(deref(rhs.thisptr))

    # def get_unit(self, UnitSpecies usp):
    #     return self.thisptr.get_unit(deref(usp.thisptr))

    def num_units(self):
        return self.thisptr.num_units()

    def deserialize(self, string serial):
        self.thisptr.deserialize(serial)

cdef Species Species_from_Cpp_Species(Cpp_Species *sp):
    cdef Cpp_Species *new_obj = new Cpp_Species(deref(sp))
    r = Species()
    del r.thisptr
    r.thisptr = new_obj
    return r

def spmatch(Species pttrn, Species sp):
    return context.spmatch(deref(pttrn.thisptr), deref(sp.thisptr))

def count_spmatches(Species pttrn, Species sp):
    return context.count_spmatches(deref(pttrn.thisptr), deref(sp.thisptr))

def rrmatch(ReactionRule pttrn, reactants):
    cdef vector[Cpp_Species] cpp_reactants
    for sp in reactants:
        cpp_reactants.push_back(deref((<Species> sp).thisptr))
    return context.rrmatch(deref(pttrn.thisptr), cpp_reactants)

def count_rrmatches(ReactionRule pttrn, reactants):
    cdef vector[Cpp_Species] cpp_reactants
    for sp in reactants:
        cpp_reactants.push_back(deref((<Species> sp).thisptr))
    return context.count_rrmatches(deref(pttrn.thisptr), cpp_reactants)

def rrgenerate(ReactionRule pttrn, reactants):
    cdef vector[Cpp_Species] cpp_reactants
    for sp in reactants:
        cpp_reactants.push_back(deref((<Species> sp).thisptr))
    cdef vector[vector[Cpp_Species]] cpp_products_list = \
        context.rrgenerate(deref(pttrn.thisptr), cpp_reactants)
    cdef vector[vector[Cpp_Species]].iterator it1 = cpp_products_list.begin()
    cdef vector[Cpp_Species].iterator it2
    retval = []
    while it1 != cpp_products_list.end():
        retval.append([])
        it2 = deref(it1).begin()
        while it2 != deref(it1).end():
            retval[-1].append(Species_from_Cpp_Species(address(deref(it2))))
            inc(it2)
        inc(it1)
    return retval

