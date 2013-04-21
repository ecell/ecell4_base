from cython.operator cimport dereference as deref
from cython cimport address, declare
from libcpp.vector cimport vector
from libcpp.string cimport string

cdef extern from "ecell4/core/Position3.hpp" namespace "ecell4":
    cdef cppclass Position3:
        Position3() except +
        Position3(Real, Real, Real) except +
        Position3(Position3 &pos3) except +

    Position3 add(Position3, Position3)
    Position3 subtract(Position3, Position3)
    Position3 divide(Position3, Real)
    Position3 multiply(Position3, Real)
    Position3 modulo(Position3, Real)
    Position3 modulo(Position3, Position3)
    Position3 abs(Position3)
    Real dot_product(Position3, Position3)
    Position3 cross_product(Position3, Position3)
    Integer length_sq(Position3)
    Integer length(Position3)
    Position3 operator+(Position3, Position3)
    Position3 operator-(Position3, Position3)
    Position3 operator/(Position3, Real)
    Position3 operator*(Position3, Real)

cdef class PyPosition3:
    #cdef Position3 *thisptr
    def __cinit__(self, Real p1, Real p2, Real p3):
        self.thisptr = new Position3(p1, p2, p3)
    def __dealloc__(self):
        del self.thisptr

cdef to_PyPosition3(Position3 *p):
    cdef Position3 *new_obj = new Position3( <Position3> deref(p) )
    r = PyPosition3(0.0, 0.0, 0.0)
    del r.thisptr
    r.thisptr = new_obj
    return r

def PyPosition3_add(PyPosition3 p1, PyPosition3 p2):
    cdef Position3 r = add(deref(p1.thisptr), deref(p2.thisptr)) 
    return to_PyPosition3( address(r) )

def PyPosition3_subtract(PyPosition3 p1, PyPosition3 p2):
    cdef Position3 r = subtract( deref(p1.thisptr), deref(p2.thisptr) )
    return to_PyPosition3( address(r) )

def PyPosotion3_divide(PyPosition3 p1, Real p2):
    cdef Position3 r = divide( deref(p1.thisptr), p2 )
    return to_PyPosition3( address(r) )

def PyPosition3_multiply(PyPosition3 p1, Real p2):
    cdef Position3 r = multiply(deref(p1.thisptr), p2)
    return to_PyPosition3( address(r) )
'''
def PyPosition3_modulo(PyPosition3 p1, PyPosition3 p2):
    cdef Position3 r = modulo( deref(p1.thisptr), <Position3>deref(p2.thisptr) )
    return to_PyPosition3( address(r) )
'''

def PyPosition3_abs(PyPosition3 p1):
    cdef Position3 r = abs( deref(p1.thisptr) )
    return to_PyPosition3( address(r) )

def PyPosition3_dot_product(PyPosition3 p1, PyPosition3 p2):
    return dot_product( deref(p1.thisptr), deref(p2.thisptr) )

def PyPosition3_cross_product(PyPosition3 p1, PyPosition3 p2):
    cdef Position3 r = cross_product( deref(p1.thisptr), deref(p2.thisptr) )
    return to_PyPosition3( address(r) )

def PyPosition3_length_sq(PyPosition3 p1):
    return length_sq( deref(p1.thisptr) )

def PyPosition3_length(PyPosition3 p1):
    return length( deref(p1.thisptr) )
