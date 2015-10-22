from types cimport Real
from core cimport Cpp_Real3


## Cpp_Real3
#  ecell4::Real3
cdef extern from "ecell4/core/Real3.hpp" namespace "ecell4":
    Cpp_Real3 add(Cpp_Real3, Cpp_Real3)
    Cpp_Real3 subtract(Cpp_Real3, Cpp_Real3)
    Cpp_Real3 divide(Cpp_Real3, Real)
    Cpp_Real3 multiply(Cpp_Real3, Real)
    Cpp_Real3 modulo(Cpp_Real3, Real)
    Cpp_Real3 modulo(Cpp_Real3, Cpp_Real3)
    Cpp_Real3 abs(Cpp_Real3)
    Real dot_product(Cpp_Real3, Cpp_Real3)
    Cpp_Real3 cross_product(Cpp_Real3, Cpp_Real3)
    Real length_sq(Cpp_Real3)
    Real length(Cpp_Real3)
