from cython.operator cimport dereference as deref, preincrement as inc
from cython cimport address

cimport create_reaction_rule as crr


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
        cdef vector[Cpp_Species] reactants = self.thisptr.reactants()
        retval = []
        cdef vector[Cpp_Species].iterator it = reactants.begin()
        while it != reactants.end():
            retval.append(
                Species_from_Cpp_Species(<Cpp_Species*>address(deref(it))))
            inc(it)
        return retval

    def products(self):
        cdef vector[Cpp_Species] products = self.thisptr.products()
        retval = []
        cdef vector[Cpp_Species].iterator it = products.begin()
        while it != products.end():
            retval.append(
                Species_from_Cpp_Species(<Cpp_Species*>address(deref(it))))
            inc(it)
        return retval

    def add_reactant(self, Species sp):
        self.thisptr.add_reactant(deref(sp.thisptr))

    def add_product(self, Species sp):
        self.thisptr.add_product(deref(sp.thisptr))

    def as_string(self):
        return self.thisptr.as_string().decode('UTF-8')

    def count(self, reactants):
        cdef vector[Cpp_Species] cpp_reactants
        for sp in reactants:
            cpp_reactants.push_back(deref((<Species> sp).thisptr))
        return self.thisptr.count(cpp_reactants)

    def set_ratelaw(self, ratelaw):
        if (isinstance(ratelaw, RatelawMassAction)):
            self.set_ratelaw_massaction(ratelaw)
        else:
            pass

    def set_ratelaw_massaction(self, RatelawMassAction ratelaw):
        self.thisptr.set_ratelaw(deref(ratelaw.thisptr))

    def generate(self, reactants):
        cdef vector[Cpp_Species] cpp_reactants
        for sp in reactants:
            cpp_reactants.push_back(deref((<Species> sp).thisptr))
        cdef vector[Cpp_ReactionRule] cpp_rules = self.thisptr.generate(cpp_reactants)
        cdef vector[Cpp_ReactionRule].iterator it1 = cpp_rules.begin()
        retval = []
        while it1 != cpp_rules.end():
            retval.append(ReactionRule_from_Cpp_ReactionRule(address(deref(it1))))
            inc(it1)
        return retval

cdef ReactionRule ReactionRule_from_Cpp_ReactionRule(Cpp_ReactionRule *rr):
    cdef Cpp_ReactionRule *new_obj = new Cpp_ReactionRule(deref(rr))
    r = ReactionRule()
    del r.thisptr
    r.thisptr = new_obj
    return r

def create_degradation_reaction_rule(Species reactant1, Real k):
    cdef Cpp_ReactionRule rr = crr.create_degradation_reaction_rule(
        deref(reactant1.thisptr), k)
    return ReactionRule_from_Cpp_ReactionRule(address(rr))

def create_synthesis_reaction_rule(Species product1, Real k):
    cdef Cpp_ReactionRule rr = crr.create_synthesis_reaction_rule(
        deref(product1.thisptr), k)
    return ReactionRule_from_Cpp_ReactionRule(address(rr))

def create_unimolecular_reaction_rule(Species reactant1, Species product1, Real k):
    cdef Cpp_ReactionRule rr = crr.create_unimolecular_reaction_rule(
        deref(reactant1.thisptr), deref(product1.thisptr), k)
    return ReactionRule_from_Cpp_ReactionRule(address(rr))

def create_binding_reaction_rule(
    Species reactant1, Species reactant2, Species product1, Real k):
    cdef Cpp_ReactionRule rr = crr.create_binding_reaction_rule(
        deref(reactant1.thisptr), deref(reactant2.thisptr),
        deref(product1.thisptr), k)
    return ReactionRule_from_Cpp_ReactionRule(address(rr))

def create_unbinding_reaction_rule(
    Species reactant1, Species product1, Species product2, Real k):
    cdef Cpp_ReactionRule rr = crr.create_unbinding_reaction_rule(
        deref(reactant1.thisptr),
        deref(product1.thisptr), deref(product2.thisptr), k)
    return ReactionRule_from_Cpp_ReactionRule(address(rr))

# def create_repulsive_reaction_rule(Species reactant1, Species reactant2):
#     cdef Cpp_ReactionRule rr = crr.create_repulsive_reaction_rule(
#         deref(reactant1.thisptr), deref(reactant2.thisptr))
#     return ReactionRule_from_Cpp_ReactionRule(address(rr))
