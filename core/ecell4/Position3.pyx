from cython.operator cimport dereference as deref
from cython cimport address


cdef class Position3:

    def __cinit__(self, Real p1, Real p2, Real p3):
        self.thisptr = new Cpp_Position3(p1, p2, p3)

    def __dealloc__(self):
        del self.thisptr

cdef Position3 Cpp_Position3_to_Position3(Cpp_Position3 *p):
    cdef Cpp_Position3 *new_obj = new Cpp_Position3(<Cpp_Position3> deref(p))
    r = Position3(0.0, 0.0, 0.0)
    del r.thisptr
    r.thisptr = new_obj
    return r

def Position3_add(Position3 p1, Position3 p2):
    cdef Cpp_Position3 r = add(deref(p1.thisptr), deref(p2.thisptr))
    return Cpp_Position3_to_Position3(address(r))

def Position3_subtract(Position3 p1, Position3 p2):
    cdef Cpp_Position3 r = subtract(deref(p1.thisptr), deref(p2.thisptr))
    return Cpp_Position3_to_Position3(address(r))

def PyPosotion3_divide(Position3 p1, Real p2):
    cdef Cpp_Position3 r = divide(deref(p1.thisptr), p2)
    return Cpp_Position3_to_Position3(address(r))

def Position3_multiply(Position3 p1, Real p2):
    cdef Cpp_Position3 r = multiply(deref(p1.thisptr), p2)
    return Cpp_Position3_to_Position3(address(r))

'''
def Position3_modulo(Position3 p1, Position3 p2):
    cdef Cpp_Position3 r = modulo(deref(p1.thisptr), <Position3>deref(p2.thisptr))
    return Cpp_Position3_to_Position3(address(r))
'''

def Position3_abs(Position3 p1):
    cdef Cpp_Position3 r = abs(deref(p1.thisptr))
    return Cpp_Position3_to_Position3(address(r))

def Position3_dot_product(Position3 p1, Position3 p2):
    return dot_product(deref(p1.thisptr), deref(p2.thisptr))

def Position3_cross_product(Position3 p1, Position3 p2):
    cdef Cpp_Position3 r = cross_product(deref(p1.thisptr), deref(p2.thisptr))
    return Cpp_Position3_to_Position3(address(r))

def Position3_length_sq(Position3 p1):
    return length_sq(deref(p1.thisptr))

def Position3_length(Position3 p1):
    return length(deref(p1.thisptr))
