from cython.operator cimport dereference as deref, preincrement as inc
from libcpp.string cimport string
from cython cimport address
cimport util

cimport context


cdef class Species:

    def __cinit__(self, serial=None, radius=None, D=None, location=None):
        if serial is None:
            self.thisptr = new Cpp_Species()
        elif radius is not None and D is not None:
            if location is None:
                self.thisptr = new Cpp_Species(
                    tostring(serial),
                    tostring(radius),
                    tostring(D))
            else:
                self.thisptr = new Cpp_Species(
                    tostring(serial),
                    tostring(radius),
                    tostring(D),
                    tostring(location))
        else:
            self.thisptr = new Cpp_Species(tostring(serial)) #XXX:

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
        return hash(self.thisptr.serial().decode('UTF-8'))

    def serial(self):
        return self.thisptr.serial().decode('UTF-8')

    def get_attribute(self, attr_name):
        return self.thisptr.get_attribute(
            tostring(attr_name)).decode('UTF-8')

    def set_attribute(self, name, value):
        self.thisptr.set_attribute(tostring(name), tostring(value))

    def remove_attribute(self, name):
        self.thisptr.remove_attribute(tostring(name))

    def has_attribute(self, name):
        return self.thisptr.has_attribute(tostring(name))

    def list_attributes(self):
        retval = self.thisptr.list_attributes()
        return [(key.decode('UTF-8'), value.decode('UTF-8'))
            for key, value in retval]

    def add_unit(self, UnitSpecies usp):
        self.thisptr.add_unit(deref(usp.thisptr))

    def count(self, Species pttrn):
        return self.thisptr.count(deref(pttrn.thisptr))

    # def get_unit(self, UnitSpecies usp):
    #     return self.thisptr.get_unit(deref(usp.thisptr))

    def units(self):
        cdef vector[Cpp_UnitSpecies] usps = self.thisptr.units()
        retval = []
        cdef vector[Cpp_UnitSpecies].iterator it = usps.begin()
        while it != usps.end():
            retval.append(UnitSpecies_from_Cpp_UnitSpecies(
            <Cpp_UnitSpecies*>(address(deref(it)))))
            inc(it)
        return retval

    def num_units(self):
        return self.thisptr.num_units()

    def deserialize(self, serial):
        self.thisptr.deserialize(tostring(serial))

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
