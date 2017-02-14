cdef class Space:
    """An abstract base class of all worlds. This is for developers.

    Space()

    """

    def __init__(self):
        """Constructor"""
        pass

    def __cinit__(self):
        self.thisptr = new shared_ptr[Cpp_Space](
            <Cpp_Space*>(new Cpp_CompartmentSpaceVectorImpl(
                Cpp_Real3(1, 1, 1)))) #XXX: DUMMY

    def __dealloc__(self):
        del self.thisptr

    def volume(self):
        return self.thisptr.get().volume()
    
    def t(self):
        return self.thisptr.get().t()

    def num_molecules(self, Species sp):
        return self.thisptr.get().num_molecules(deref(sp.thisptr))


cdef Space Space_from_Cpp_Space(shared_ptr[Cpp_Space] space):
    r = Space()
    r.thisptr.swap(space)
    return r
