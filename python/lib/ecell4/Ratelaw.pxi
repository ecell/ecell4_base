
from cython.operator cimport dereference as deref, preincrement as inc
from cython cimport address

from libcpp.vector cimport vector

cdef class RatelawMassAction:
    def __cinit__(self, Real k):
        self.thisptr = new shared_ptr[Cpp_RatelawMassAction]( <Cpp_RatelawMassAction*>(new Cpp_RatelawMassAction(k)))
    def __dealloc__(self):
        del self.thisptr
    def set_k(self, Real k):
        self.thisptr.get().set_k(k)
    def get_k(self):
        return self.thisptr.get().get_k()


cdef double indirect_function(void *func, vector[Real] reactants, vector[Real] products, Real volume):
    py_reactants = []
    cdef vector[Real].iterator it1 = reactants.begin()
    while it1 != reactants.end():
        py_reactants.append(deref(it1))
        inc(it1)
    py_products = []
    cdef vector[Real].iterator it2 = products.begin()
    while it2 != products.end():
        py_products.append(deref(it1))
        inc(it2)
    return (<object>func)(py_reactants, py_products, volume)

cdef class RatelawCallback:
    def __cinit__(self, pyfunc):
        self.thisptr = new shared_ptr[Cpp_RatelawCythonCallback]( <Cpp_RatelawCythonCallback*>(new Cpp_RatelawCythonCallback(<Indirect_Functype>indirect_function, <void*>pyfunc)) )

    def __dealloc__(self):
        del self.thisptr
    def set_callback(self, pyfunc):
        self.thisptr.get().set_callback_pyfunc(<Python_Functype>pyfunc)

    def call(self):
        return self.thisptr.get().call()
