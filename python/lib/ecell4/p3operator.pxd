from types cimport Real
from core cimport Cpp_Position3


## Cpp_Position3
#  ecell4::Position3
cdef extern from "ecell4/core/Position3.hpp" namespace "ecell4":
    Cpp_Position3 add(Cpp_Position3, Cpp_Position3)
    Cpp_Position3 subtract(Cpp_Position3, Cpp_Position3)
    Cpp_Position3 divide(Cpp_Position3, Real)
    Cpp_Position3 multiply(Cpp_Position3, Real)
    Cpp_Position3 modulo(Cpp_Position3, Real)
    Cpp_Position3 modulo(Cpp_Position3, Cpp_Position3)
    Cpp_Position3 abs(Cpp_Position3)
    Real dot_product(Cpp_Position3, Cpp_Position3)
    Cpp_Position3 cross_product(Cpp_Position3, Cpp_Position3)
    Real length_sq(Cpp_Position3)
    Real length(Cpp_Position3)
