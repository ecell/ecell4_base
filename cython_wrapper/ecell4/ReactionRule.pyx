from cython.operator cimport dereference as deref, preincrement as inc
from cython cimport address


cdef class ReactionRule:

    def __cinit__(self):
        self.thisptr = new Cpp_ReactionRule()

    def __dealloc__(self):
        del self.thisptr

    def k(self):
        return self.thisptr.k()

    def set_k(self, Real k):
        self.thisptr.set_k(k)

    def reactants(self):
        cdef multiset[Cpp_Species] reactants = self.thisptr.reactants()
        retval = []
        cdef multiset[Cpp_Species].iterator it = reactants.begin()
        while it != reactants.end():
            retval.append(
                Cpp_Species_to_Species(<Cpp_Species*>address(deref(it))))
            inc(it)
        return retval

    def products(self):
        cdef multiset[Cpp_Species] products = self.thisptr.products()
        retval = []
        cdef multiset[Cpp_Species].iterator it = products.begin()
        while it != products.end():
            retval.append(
                Cpp_Species_to_Species(<Cpp_Species*>address(deref(it))))
            inc(it)
        return retval

    def add_reactant(self, Species sp):
        self.thisptr.add_reactant(deref(sp.thisptr))

    def add_product(self, Species sp):
        self.thisptr.add_product(deref(sp.thisptr))

def create_unimolecular_reaction_rule(
        Species reactant1, Species product1, Real k):
    rr = ReactionRule()
    rr.add_reactant(reactant1)
    rr.add_product(product1)
    rr.set_k(k)
    return rr

def create_binding_reaction_rule(
        Species reactant1, Species reactant2, Species product1, Real k):
    rr = ReactionRule()
    rr.add_reactant(reactant1)
    rr.add_reactant(reactant2)
    rr.add_product(product1)
    rr.set_k(k)
    return rr

def create_unbinding_reaction_rue(
        Species reactant1, Species product1, Species product2, Real k):
    rr = ReactionRule()
    rr.add_reactant(reactant1)
    rr.add_product(product1)
    rr.add_product(product2)
    rr.set_k(k)
    return rr

def create_repulsive_reaction_rule(
        Species reactant1, Species reactant2):
    rr = ReactionRule()
    rr.set_k(0.0)
    rr.add_reactant(reactant1)
    rr.add_product(reactant2)
    return rr
