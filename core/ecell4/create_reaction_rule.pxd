from types cimport Real
from core cimport Cpp_Species, Cpp_ReactionRule


## utility functions in Model.hpp
#  ecell4::Model
cdef extern from "ecell4/core/Model.hpp" namespace "ecell4":
    Cpp_ReactionRule create_degradation_reaction_rule(Cpp_Species&, Real)
    Cpp_ReactionRule create_synthesis_reaction_rule(Cpp_Species&, Real)
    Cpp_ReactionRule create_unimolecular_reaction_rule(
        Cpp_Species&, Cpp_Species&, Real)
    Cpp_ReactionRule create_binding_reaction_rule(
        Cpp_Species&, Cpp_Species&, Cpp_Species&, Real)
    Cpp_ReactionRule create_unbinding_reaction_rule(
        Cpp_Species&, Cpp_Species&, Cpp_Species&, Real)
    Cpp_ReactionRule create_repulsive_reaction_rule(
        Cpp_Species&, Cpp_Species&)

