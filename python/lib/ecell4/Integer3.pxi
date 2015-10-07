from cython.operator cimport dereference as deref
from cython cimport address


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

cdef Integer3 Integer3_from_Cpp_Integer3(Cpp_Integer3 *p):
    cdef Cpp_Integer3 *new_obj = new Cpp_Integer3(<Cpp_Integer3> deref(p))
    r = Integer3(0.0, 0.0, 0.0)
    del r.thisptr
    r.thisptr = new_obj
    return r
