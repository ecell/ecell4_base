from cython.operator cimport dereference as deref
from cython cimport address
cimport p3operator


cdef class Real3:

    def __cinit__(self, Real p1, Real p2, Real p3):
        self.thisptr = new Cpp_Real3(p1, p2, p3)

    def __dealloc__(self):
        del self.thisptr

    def __getitem__(self, Integer i):
        if i > 2:
            raise IndexError("index out of bounds")
        return deref(self.thisptr)[i]

    def __add__(Real3 self, Real3 other):
        return add(self, other)

    def __sub__(Real3 self, Real3 other):
        return subtract(self, other)

    def __div__(Real3 self, Real other):
        return divide(self, other)

    def __truediv__(Real3 self, Real other):
        return divide(self, other)

    def __mul__(self, other):
        print("Real3.__mul__({0}, {1})".format(self, other))
        if isinstance(self, Real3):
            return multiply(<Real3>self, <Real>other)
        elif isinstance(other, Real3):
            return multiply(<Real3>other, <Real>self)
        else:
            raise ValueError(
                'invalid value was given: '
                + repr(self) + ' : ' + repr(other))


cdef Real3 Real3_from_Cpp_Real3(Cpp_Real3 *p):
    cdef Cpp_Real3 *new_obj = new Cpp_Real3(<Cpp_Real3> deref(p))
    r = Real3(0.0, 0.0, 0.0)
    del r.thisptr
    r.thisptr = new_obj
    return r

def add(Real3 p1, Real3 p2):
    cdef Cpp_Real3 r = p3operator.add(deref(p1.thisptr), deref(p2.thisptr))
    return Real3_from_Cpp_Real3(address(r))

def subtract(Real3 p1, Real3 p2):
    cdef Cpp_Real3 r = p3operator.subtract(deref(p1.thisptr), deref(p2.thisptr))
    return Real3_from_Cpp_Real3(address(r))

def divide(Real3 p1, Real p2):
    cdef Cpp_Real3 r = p3operator.divide(deref(p1.thisptr), p2)
    return Real3_from_Cpp_Real3(address(r))

def multiply(Real3 p1, Real p2):
    cdef Cpp_Real3 r = p3operator.multiply(deref(p1.thisptr), p2)
    return Real3_from_Cpp_Real3(address(r))

# def modulo(Real3 p1, Real3 p2):
#     cdef Cpp_Real3 r = p3operator.modulo(
#         deref(p1.thisptr), <Real3>deref(p2.thisptr))
#     return Real3_from_Cpp_Real3(address(r))

def abs(Real3 p1):
    cdef Cpp_Real3 r = p3operator.abs(deref(p1.thisptr))
    return Real3_from_Cpp_Real3(address(r))

def dot_product(Real3 p1, Real3 p2):
    return p3operator.dot_product(deref(p1.thisptr), deref(p2.thisptr))

def cross_product(Real3 p1, Real3 p2):
    cdef Cpp_Real3 r = p3operator.cross_product(deref(p1.thisptr), deref(p2.thisptr))
    return Real3_from_Cpp_Real3(address(r))

def length_sq(Real3 p1):
    return p3operator.length_sq(deref(p1.thisptr))

def length(Real3 p1):
    return p3operator.length(deref(p1.thisptr))
