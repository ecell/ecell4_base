cdef class GSLRandomNumberGenerator:

    def __cinit__(self):
        # self.thisptr = new shared_ptr[Cpp_GSLRandomNumberGenerator](
        #     new Cpp_GSLRandomNumberGenerator())
        self.thisptr = new shared_ptr[Cpp_RandomNumberGenerator](
            <Cpp_RandomNumberGenerator*> (new Cpp_GSLRandomNumberGenerator()))

    def __dealloc__(self):
        del self.thisptr

    def uniform(self, Real min, Real max):
        return self.thisptr.get().uniform(min, max)

    def uniform_int(self, Integer min, Integer max):
        return self.thisptr.get().uniform_int(min, max)

    def gaussian(self, Real mean, Real sigma):
        return self.thisptr.get().gaussian(mean, sigma)

    def binomial(self, Real p, Integer n):
        return self.thisptr.get().binomial(p, n)

    def seed(self, val = None):
        if val is None:
            self.thisptr.get().seed()
        else:
            self.thisptr.get().seed(<Integer> val)

cdef GSLRandomNumberGenerator GSLRandomNumberGenerator_from_Cpp_RandomNumberGenerator(
    shared_ptr[Cpp_RandomNumberGenerator] rng):
    r = GSLRandomNumberGenerator()
    r.thisptr.swap(rng)
    return r
