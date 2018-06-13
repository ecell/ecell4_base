cdef class WorldInterface:
    """An abstract base class of all worlds. This is for developers.

    WorldInterface()

    """

    def __init__(self):
        """Constructor"""
        pass

    def __cinit__(self):
        self.thisptr = new shared_ptr[Cpp_WorldInterface](
            <Cpp_WorldInterface*>(new Cpp_CompartmentSpaceVectorImpl(
                Cpp_Real3(1, 1, 1)))) #XXX: DUMMY

    def __dealloc__(self):
        del self.thisptr

    def volume(self):
        return self.thisptr.get().volume()

    def t(self):
        return self.thisptr.get().t()

    def num_molecules(self, Species sp):
        return self.thisptr.get().num_molecules(deref(sp.thisptr))


cdef WorldInterface WorldInterface_from_Cpp_WorldInterface(shared_ptr[Cpp_WorldInterface] world):
    r = WorldInterface()
    r.thisptr.swap(world)
    return r
