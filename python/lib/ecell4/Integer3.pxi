from cython.operator cimport dereference as deref
from cython cimport address


cdef class Integer3:

    def __cinit__(self, Integer col, Integer row, Integer layer):
        self.thisptr = new Cpp_Integer3(col, row, layer)

    def __dealloc__(self):
        del self.thisptr

    @property
    def col(self):
        return self.thisptr.col

    @property
    def row(self):
        return self.thisptr.row

    @property
    def layer(self):
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
