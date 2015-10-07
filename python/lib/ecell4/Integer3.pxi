from cython.operator cimport dereference as deref
from cython cimport address
cimport integer3operators

cdef class Integer3:
    """A class representing a vector consisting of three integers.

    Integer3(Integer p1, Integer p2, Integer p3)

    """

    def __init__(self, Integer p1, Integer p2, Integer p3):
        """Constructor.

        Args:
            p1 (Integer): The first value in the vector.
            p2 (Integer): The second value in the vector.
            p3 (Integer): The third value in the vector.

        """
        pass

    def __cinit__(self, Integer col, Integer row, Integer layer):
        self.thisptr = new Cpp_Integer3(col, row, layer)

    def __dealloc__(self):
        del self.thisptr

    @property
    def col(self):
        """Return the first value."""
        return self.thisptr.col

    @property
    def row(self):
        """Return the second value."""
        return self.thisptr.row

    @property
    def layer(self):
        """Return the third value."""
        return self.thisptr.layer

    def __getitem__(self, Integer i):
        if i > 2:
            raise IndexError("index out of bounds")
        return deref(self.thisptr)[i]

    def __add__(Integer3 self, Integer3 other):
        return integer3_add(self, other)

    def __sub__(Integer3 self, Integer3 other):
        return integer3_subtract(self, other)

    def __abs__(Integer3 self):
        return integer3_abs(self)

    # def __div__(Integer3 self, Integer other):
    #     return integer3_divide(self, other)

    # def __truediv__(Integer3 self, Integer other):
    #     return integer3_divide(self, other)

    # def __mul__(self, other):
    #     if isinstance(self, Integer3):
    #         return integer3_multiply(<Integer3>self, <Integer>other)
    #     elif isinstance(other, Integer3):
    #         return integer3_multiply(<Integer3>other, <Integer>self)
    #     else:
    #         raise ValueError(
    #             'invalid value was given: '
    #             + repr(self) + ' : ' + repr(other))

cdef Integer3 Integer3_from_Cpp_Integer3(Cpp_Integer3 *p):
    cdef Cpp_Integer3 *new_obj = new Cpp_Integer3(<Cpp_Integer3> deref(p))
    r = Integer3(0.0, 0.0, 0.0)
    del r.thisptr
    r.thisptr = new_obj
    return r

def integer3_add(Integer3 p1, Integer3 p2):
    """add(p1, p2) -> Integer3

    Add two ``Integer3``s, and returns the sum.

    Args:
        p1 (Integer3): The first vector.
        p2 (Integer3): The second vector.

    Returns:
        Integer3: The sum of two vectors, ``p1 + p2``.

    """
    cdef Cpp_Integer3 r = integer3operators.add(deref(p1.thisptr), deref(p2.thisptr))
    return Integer3_from_Cpp_Integer3(address(r))

def integer3_subtract(Integer3 p1, Integer3 p2):
    """subtract(p1, p2) -> Integer3

    Subtract p2 from p1.

    Args:
        p1 (Integer3): The left-hand-side vector.
        p2 (Integer3): The right-hand-side vector.

    Returns:
        Integer3: Its difference, ``p1 - p2``.

    """
    cdef Cpp_Integer3 r = integer3operators.subtract(deref(p1.thisptr), deref(p2.thisptr))
    return Integer3_from_Cpp_Integer3(address(r))

def integer3_abs(Integer3 p1):
    """real3_abs(p1) -> Integer3

    Return an absolute vector of the given vector.

    Args:
        p1 (Integer3): A vector.

    Returns:
        Integer3: The absolute vector, which consists of absolute value
            of the given vector.

    Note:
        This is NOT for taking the norm of a vector. See ``length`` also.

    """
    cdef Cpp_Integer3 r = integer3operators.abs(deref(p1.thisptr))
    return Integer3_from_Cpp_Integer3(address(r))
