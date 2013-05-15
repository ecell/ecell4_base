from cython.operator cimport dereference as deref


cdef class NetworkModel:

    def __cinit__(self):
        # self.thisptr = new NetworkModel()
        self.thisptr = new shared_ptr[Cpp_NetworkModel](new Cpp_NetworkModel())

    def __dealloc__(self):
        del self.thisptr

    def add_species(self, Species sp):
        self.thisptr.get().add_species(deref(sp.thisptr))

    def has_species(self, Species sp):
        return self.thisptr.get().has_species(deref(sp.thisptr))

    def remove_species(self, Species sp):
        self.thisptr.get().remove_species(deref(sp.thisptr))

    def add_reaction_rule(self, ReactionRule rr):
        self.thisptr.get().add_reaction_rule(deref(rr.thisptr))

    def remove_reaction_rule(self, ReactionRule rr):
        self.thisptr.get().remove_reaction_rule(deref(rr.thisptr))

    def has_reaction_rule(self, ReactionRule rr):
        self.thisptr.get().has_reaction_rule(deref(rr.thisptr))

    # def add_reactant(self, PySpecies sp):
    #     self.thisptr.add_reactant(deref(sp.thisptr))
    # def add_product(self, PySpecies sp):
    #     self.thisptr.add_product(deref(sp.thisptr))
    # def reactants(self):
    #     # self.thisptr.reactants()
    #     pass
