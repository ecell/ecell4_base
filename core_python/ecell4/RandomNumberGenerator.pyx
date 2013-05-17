cdef class GSLRandomNumberGenerator:

    def __cinit__(self):
        self.thisptr = new shared_ptr[Cpp_GSLRandomNumberGenerator](
            new Cpp_GSLRandomNumberGenerator())

    def __dealloc__(self):
        del self.thisptr

    def uniform(self, Real min, Real max):
        return self.thisptr.get().uniform(min, max)

    def uniform_int(self, Integer min, Integer max):
        return self.thisptr.get().uniform_int(min, max)

    def gaussian(self, Real mean, Real sigma):
        return self.thisptr.get().gaussian(mean, sigma)

    def seed(self, Integer val):
        self.thisptr.get().seed(val)
