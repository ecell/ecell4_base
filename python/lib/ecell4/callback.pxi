

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


cdef bool indirect_func_space(
        void *pyfunc, shared_ptr[Cpp_CompartmentSpaceVectorImpl] sp):
    sp_obj = CompartmentSpaceVectorImpl(Real3(0., 0., 0.))
    del sp_obj.thisptr
    sp_obj.thisptr = sp.get()
    return (<object>pyfunc)(sp_obj)

cdef class PythonNumberHooker:
    def __init__(self, pyfunc):
        pass
    def __cinit__(self, pyfunc):
        self.thisptr = new Cpp_PythonNumberHooker(<stepladder_type_space>indirect_func_space, <void*>pyfunc)
        self.pyfunc_ = pyfunc
    def is_available(self):
        return self.thisptr.is_available()
    #def call(self, space):
        # This is implemented for Test and Debug.
        #return self.thisptr.call(space)

cdef class PythonCompartmentSpaceVectorImplHooker:
    def __init__(self, pyfunc):
        pass
    def __cinit__(self, pyfunc):
        self.thisptr = new PythonHook_1arg[shared_ptr[Cpp_CompartmentSpaceVectorImpl]](
                indirect_func_space, <void*>pyfunc)
        self.pyfunc_ = pyfunc

    def __dealloc__(self):
        del self.thisptr

    def is_available(self):
        return self.thisptr.is_available()

        

