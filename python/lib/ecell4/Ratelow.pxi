
from cython.operator cimport dereference as deref, preincrement as inc
from cython cimport address

from libcpp.vector cimport vector

cdef class RatelowMassAction:
    def __cinit__(self, Real k):
        self.thisptr = new shared_ptr[Cpp_RatelowMassAction]( <Cpp_RatelowMassAction*>(new Cpp_RatelowMassAction(k)))
    def __dealloc__(self):
        del self.thisptr
    def set_k(self, Real k):
        self.thisptr.get().set_k(k)
    def get_k(self):
        return self.thisptr.get().get_k()


