

cdef extern from "gsl/gsl_rng.h":
    ctypedef struct gsl_rng:
        pass
        

cdef extern from "ecell4/core/RandomNumberGenerator.hpp" namespace "ecell4":
    cdef cppclass GSLRandomNumberGenerator:
        #GSLRandomNumberGenerator(shared_ptr[gsl_rng]) except +
        GSLRandomNumberGenerator() except +
        Real uniform(Real, Real)
        Integer uniform_int(Integer, Integer)
        Real gaussian(Real, Real)
        void seed(Integer)

cdef class PyRandomNumberGenerator:
    cdef GSLRandomNumberGenerator *thisptr
    def __cinit__(self):
        self.thisptr = new GSLRandomNumberGenerator()
    def __dealloc__(self):
        del self.thisptr

    def uniform(self, Real min, Real max):
        return self.thisptr.uniform(min, max)
    def uniform_inf(self, Integer min, Integer max):
        return self.thisptr.uniform_int(min, max)
    def gaussian(self, Real mean, Real sigma):
        return self.thisptr.gaussian(mean, sigma)
    def seed(self, Integer val):
        self.thisptr.seed(val)
