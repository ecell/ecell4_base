from types cimport *


cdef extern from "ecell4/core/functions.hpp" namespace "ecell4":
    Real cbrt(const Real&)

