from cython.operator cimport dereference as deref
from cython cimport address


cdef class Global:

    def __cinit__(self, Integer col, Integer row, Integer layer):
        self.thisptr = new Cpp_Global(col, row, layer)

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

cdef Global Global_from_Cpp_Global(Cpp_Global *p):
    cdef Cpp_Global *new_obj = new Cpp_Global(<Cpp_Global> deref(p))
    r = Global(0.0, 0.0, 0.0)
    del r.thisptr
    r.thisptr = new_obj
    return r
