from libcpp cimport bool
from libcpp.vector cimport vector

from types cimport Integer
from core cimport Cpp_Species, Cpp_ReactionRule


cdef extern from "ecell4/core/Context.hpp" namespace "ecell4":
    bool spmatch(Cpp_Species, Cpp_Species)
    Integer count_spmatches(Cpp_Species, Cpp_Species)
    bool rrmatch(Cpp_ReactionRule, vector[Cpp_Species])
    Integer count_rrmatches(Cpp_ReactionRule, vector[Cpp_Species])
    vector[vector[Cpp_Species]] rrgenerate(Cpp_ReactionRule, vector[Cpp_Species])
