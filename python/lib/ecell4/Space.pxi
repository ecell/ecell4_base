cdef class Space:

    def __cinit__(self):
        self.thisptr = new shared_ptr[Cpp_Space](
            <Cpp_Space*>(new Cpp_CompartmentSpaceVectorImpl(Cpp_Position3(1, 1, 1)))) #XXX: DUMMY

    def __dealloc__(self):
        del self.thisptr
