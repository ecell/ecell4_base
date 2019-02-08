from types cimport Integer, Real
from core cimport Cpp_Integer3


## Cpp_Integer3
#  ecell4::Integer3
cdef extern from "ecell4/core/Integer3.hpp" namespace "ecell4":
    Cpp_Integer3 add(Cpp_Integer3, Cpp_Integer3)
    Cpp_Integer3 subtract(Cpp_Integer3, Cpp_Integer3)
    # Cpp_Integer3 divide(Cpp_Integer3, Integer)
    Cpp_Integer3 multiply(Cpp_Integer3, Integer)
    # Cpp_Integer3 modulo(Cpp_Integer3, Integer)
    # Cpp_Integer3 modulo(Cpp_Integer3, Cpp_Integer3)
    Cpp_Integer3 abs(Cpp_Integer3)
    Integer dot_product(Cpp_Integer3, Cpp_Integer3)
    # Cpp_Integer3 cross_product(Cpp_Integer3, Cpp_Integer3)
    Integer length_sq(Cpp_Integer3)
    Real length(Cpp_Integer3)
