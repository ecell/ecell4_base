from cython.operator cimport dereference as deref, preincrement as inc
from cython cimport address

from libcpp.vector cimport vector


cdef class NetworkModel:

    def __cinit__(self):
        # self.thisptr = new NetworkModel()
        self.thisptr = new shared_ptr[Cpp_NetworkModel](new Cpp_NetworkModel())

    def __dealloc__(self):
        del self.thisptr

    def add_species_attribute(self, Species sp):
        self.thisptr.get().add_species_attribute(deref(sp.thisptr))

    def has_species_attribute(self, Species sp):
        return self.thisptr.get().has_species_attribute(deref(sp.thisptr))

    def remove_species_attribute(self, Species sp):
        self.thisptr.get().remove_species_attribute(deref(sp.thisptr))

    def add_reaction_rule(self, ReactionRule rr):
        self.thisptr.get().add_reaction_rule(deref(rr.thisptr))

    def remove_reaction_rule(self, ReactionRule rr):
        self.thisptr.get().remove_reaction_rule(deref(rr.thisptr))

    def has_reaction_rule(self, ReactionRule rr):
        self.thisptr.get().has_reaction_rule(deref(rr.thisptr))

    def list_species(self):
        cdef vector[Cpp_Species] species = self.thisptr.get().list_species()
        retval = []
        cdef vector[Cpp_Species].iterator it = species.begin()
        while it != species.end():
            retval.append(Species_from_Cpp_Species(
                <Cpp_Species*>(address(deref(it)))))
            inc(it)
        return retval

    def query_reaction_rules(self, Species sp1, Species sp2 = None):
        cdef vector[Cpp_ReactionRule] rules
        if sp2 is None:
            rules = self.thisptr.get().query_reaction_rules(
                deref(sp1.thisptr))
        else:
            rules = self.thisptr.get().query_reaction_rules(
                deref(sp1.thisptr), deref(sp2.thisptr))
        retval = []
        cdef vector[Cpp_ReactionRule].iterator it = rules.begin()
        while it != rules.end():
            retval.append(ReactionRule_from_Cpp_ReactionRule(
                <Cpp_ReactionRule*>(address(deref(it)))))
            inc(it)
        return retval

    # def add_reactant(self, PySpecies sp):
    #     self.thisptr.add_reactant(deref(sp.thisptr))
    # def add_product(self, PySpecies sp):
    #     self.thisptr.add_product(deref(sp.thisptr))
    # def reactants(self):
    #     # self.thisptr.reactants()
    #     pass
