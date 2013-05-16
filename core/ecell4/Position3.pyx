from cython.operator cimport dereference as deref
from cython cimport address
cimport p3operator


cdef class Position3:

    def __cinit__(self, Real p1, Real p2, Real p3):
        self.thisptr = new Cpp_Position3(p1, p2, p3)

    def __dealloc__(self):
        del self.thisptr

    def __getitem__(self, Integer i):
        if i > 2:
            raise IndexError("index out of bounds")
        return deref(self.thisptr)[i]

cdef Position3 Position3_from_Cpp_Position3(Cpp_Position3 *p):
    cdef Cpp_Position3 *new_obj = new Cpp_Position3(<Cpp_Position3> deref(p))
    r = Position3(0.0, 0.0, 0.0)
    del r.thisptr
    r.thisptr = new_obj
    return r

def add(Position3 p1, Position3 p2):
    cdef Cpp_Position3 r = p3operator.add(deref(p1.thisptr), deref(p2.thisptr))
    return Position3_from_Cpp_Position3(address(r))

def subtract(Position3 p1, Position3 p2):
    cdef Cpp_Position3 r = p3operator.subtract(deref(p1.thisptr), deref(p2.thisptr))
    return Position3_from_Cpp_Position3(address(r))

def divide(Position3 p1, Real p2):
    cdef Cpp_Position3 r = p3operator.divide(deref(p1.thisptr), p2)
    return Position3_from_Cpp_Position3(address(r))

def multiply(Position3 p1, Real p2):
    cdef Cpp_Position3 r = p3operator.multiply(deref(p1.thisptr), p2)
    return Position3_from_Cpp_Position3(address(r))

# def modulo(Position3 p1, Position3 p2):
#     cdef Cpp_Position3 r = p3operator.modulo(
#         deref(p1.thisptr), <Position3>deref(p2.thisptr))
#     return Position3_from_Cpp_Position3(address(r))

def abs(Position3 p1):
    cdef Cpp_Position3 r = p3operator.abs(deref(p1.thisptr))
    return Position3_from_Cpp_Position3(address(r))

def dot_product(Position3 p1, Position3 p2):
    return p3operator.dot_product(deref(p1.thisptr), deref(p2.thisptr))

def cross_product(Position3 p1, Position3 p2):
    cdef Cpp_Position3 r = p3operator.cross_product(deref(p1.thisptr), deref(p2.thisptr))
    return Position3_from_Cpp_Position3(address(r))

def length_sq(Position3 p1):
    return p3operator.length_sq(deref(p1.thisptr))

def length(Position3 p1):
    return p3operator.length(deref(p1.thisptr))
