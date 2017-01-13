
#from cython.operator cimport dereference as deref, preincrement as inc
#from cython cimport address
#
#from libcpp.vector cimport vector

cdef int indirect_function(
        void *pyfunc, vector[Real] aaa):

    v = []
    cdef vector[Real].iterator it1 = aaa.begin()
    while it1 != aaa.end():
        v.append(deref(it1))
        inc(it1)

    (<object>pyfunc)(v)
    return 0
    
cdef class CallbackWrapper:
    def __init__(self, pyfunc):
        pass
    def __cinit__(self, pyfunc):
        self.thisptr = new Cpp_CallbackWrapper(<stepladder_type>indirect_function, <void*>pyfunc)
        self.pyfunc_ = pyfunc

    def __dealloc__(self):
        del self.thisptr

    def is_available(self):
        """Check if this ratelaw is available or not. Return True always."""
        return self.thisptr.is_available()

    def call(self):
        self.thisptr.call()



