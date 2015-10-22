from cython.operator cimport dereference as deref
from cython cimport address
cimport real3operators


cdef class Real3:
    """A class representing a three-dimensional vector or position.

    Real3(Real p1, Real p2, Real p3)

    """

    def __init__(self, Real p1, Real p2, Real p3):
        """Constructor.

        Parameters
        ----------
        p1 : Real
            The first value in the vector.
        p2 : Real
            The second value in the vector.
        p3 : Real
            The third value in the vector.

        """
        pass

    def __cinit__(self, Real p1, Real p2, Real p3):
        self.thisptr = new Cpp_Real3(p1, p2, p3)

    def __dealloc__(self):
        del self.thisptr

    def __getitem__(self, Integer i):
        if i > 2:
            raise IndexError("index out of bounds")
        return deref(self.thisptr)[i]

    def __add__(Real3 self, Real3 other):
        return real3_add(self, other)

    def __sub__(Real3 self, Real3 other):
        return real3_subtract(self, other)

    def __div__(Real3 self, Real other):
        return real3_divide(self, other)

    def __truediv__(Real3 self, Real other):
        return real3_divide(self, other)

    def __mul__(self, other):
        if isinstance(self, Real3):
            return real3_multiply(<Real3>self, <Real>other)
        elif isinstance(other, Real3):
            return real3_multiply(<Real3>other, <Real>self)
        else:
            raise ValueError(
                'invalid value was given: '
                + repr(self) + ' : ' + repr(other))

    def __abs__(self):
        return real3_abs(self)


cdef Real3 Real3_from_Cpp_Real3(Cpp_Real3 *p):
    cdef Cpp_Real3 *new_obj = new Cpp_Real3(<Cpp_Real3> deref(p))
    r = Real3(0.0, 0.0, 0.0)
    del r.thisptr
    r.thisptr = new_obj
    return r

def real3_add(Real3 p1, Real3 p2):
    """real3_add(p1, p2) -> Real3

    Add two ``Real3``s, and returns the sum.

    Parameters
    ----------
    p1 : Real3
        The first vector.
    p2 : Real3
        The second vector.

    Returns
    -------
    Real3:
        The sum of two vectors, ``p1 + p2``.

    """
    cdef Cpp_Real3 r = real3operators.add(deref(p1.thisptr), deref(p2.thisptr))
    return Real3_from_Cpp_Real3(address(r))

def real3_subtract(Real3 p1, Real3 p2):
    """real3_subtract(p1, p2) -> Real3

    Subtract p2 from p1.

    Parameters
    ----------
    p1 : Real3
        The left-hand-side vector.
    p2 : Real3
        The right-hand-side vector.

    Returns
    -------
    Real3:
        Its difference, ``p1 - p2``.

    """
    cdef Cpp_Real3 r = real3operators.subtract(deref(p1.thisptr), deref(p2.thisptr))
    return Real3_from_Cpp_Real3(address(r))

def real3_divide(Real3 p1, Real p2):
    """real3_divide(p1, p2) -> Real3

    Divide p1 by p2.

    Parameters
    ----------
    p1 : Real3
        The numerator.
    p2 : Real
        The denominator.

    Returns
    -------
    Real3:
        The divided vector, ``p1 / p2``.

    """
    cdef Cpp_Real3 r = real3operators.divide(deref(p1.thisptr), p2)
    return Real3_from_Cpp_Real3(address(r))

def real3_multiply(Real3 p1, Real p2):
    """real3_multiply(p1, p2) -> Real3

    Multiply p1 by p2.

    Parameters
    ----------
    p1 : Real3
        A vector.
    p2 : Real
        A factor.

    Returns
    -------
    Real3:
        The multipled vector, ``p1 * p2``.

    """
    cdef Cpp_Real3 r = real3operators.multiply(deref(p1.thisptr), p2)
    return Real3_from_Cpp_Real3(address(r))

# def real3_modulo(Real3 p1, Real3 p2):
#     cdef Cpp_Real3 r = real3operators.modulo(
#         deref(p1.thisptr), <Real3>deref(p2.thisptr))
#     return Real3_from_Cpp_Real3(address(r))

def real3_abs(Real3 p1):
    """real3_abs(p1) -> Real3

    Return an absolute vector of the given vector.

    Parameters
    ----------
    p1 : Real3
        A vector.

    Returns
    -------
    Real3:
        The absolute vector, which consists of absolute value of the given vector.

    Notes
    -----
    This is NOT for taking the norm of a vector.

    See Also
    --------
    length

    """
    cdef Cpp_Real3 r = real3operators.abs(deref(p1.thisptr))
    return Real3_from_Cpp_Real3(address(r))

def dot_product(Real3 p1, Real3 p2):
    """dot_product(p1, p2) -> Real

    Return a dot product between two vectors

    """
    return real3operators.dot_product(deref(p1.thisptr), deref(p2.thisptr))

def cross_product(Real3 p1, Real3 p2):
    """cross_product(p1, p2) -> Real3

    Return a cross product between two vectors

    """
    cdef Cpp_Real3 r = real3operators.cross_product(deref(p1.thisptr), deref(p2.thisptr))
    return Real3_from_Cpp_Real3(address(r))

def length_sq(Real3 p1):
    """length_sq(p1) -> Real

    Return a square of a Euclidean norm of the given vector.

    """
    return real3operators.length_sq(deref(p1.thisptr))

def length(Real3 p1):
    """length(p1) -> Real

    Return a Euclidean norm of the given vector.
    This is almost equivalent to call ``sqrt(length_sq(p1))``

    """
    return real3operators.length(deref(p1.thisptr))
