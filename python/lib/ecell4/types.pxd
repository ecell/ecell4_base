from libc.stdint cimport int64_t

ctypedef int64_t Integer
ctypedef double Real

cdef extern from "ecell4/core/types.hpp" namespace "ecell4":
    # cdef Real inf
    cdef Real epsilon
    cdef Real N_A
