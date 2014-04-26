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
        return self.thisptr.get().has_reaction_rule(deref(rr.thisptr))

    def num_reaction_rules(self):
        return self.thisptr.get().num_reaction_rules()

    def apply_species_attributes(self, Species sp):
        cdef Cpp_Species retval = self.thisptr.get().apply_species_attributes(
            deref(sp.thisptr))
        return Species_from_Cpp_Species(address(retval))

    def create_species(self, string name):
        cdef Cpp_Species retval = self.thisptr.get().create_species(name)
        return Species_from_Cpp_Species(address(retval))

    def reaction_rules(self):
        cdef vector[Cpp_ReactionRule] c_rr_vector = self.thisptr.get().reaction_rules()
        retval = []
        cdef vector[Cpp_ReactionRule].iterator it = c_rr_vector.begin()
        while it != c_rr_vector.end():
            retval.append(ReactionRule_from_Cpp_ReactionRule(
                <Cpp_ReactionRule*>(address(deref(it)))) )
            inc(it)
        return retval

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
