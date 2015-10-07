cdef class GSLRandomNumberGenerator:
    """A random number generator using the GNU Scientific Library (GSL).

    GSLRandomNumberGenerator()

    """

    def __init__(self):
        """Constructor."""
        pass

    def __cinit__(self):
        # self.thisptr = new shared_ptr[Cpp_GSLRandomNumberGenerator](
        #     new Cpp_GSLRandomNumberGenerator())
        self.thisptr = new shared_ptr[Cpp_RandomNumberGenerator](
            <Cpp_RandomNumberGenerator*> (new Cpp_GSLRandomNumberGenerator()))

    def __dealloc__(self):
        del self.thisptr

    def uniform(self, Real min, Real max):
        """uniform(min, max) -> Real

        Return a uniform random number within the given range.

        Args:
            min (Real): The minimum value in the range.
            max (Real): The maximum value in the range.

        Returns:
            Real: A random number uniformly distributed in the range [min, max).

        """
        return self.thisptr.get().uniform(min, max)

    def uniform_int(self, Integer min, Integer max):
        """uniform_int(min, max) -> Integer

        Return a uniform random number within the given range.

        Args:
            min (Real): The minimum value in the range.
            max (Real): The maximum value in the range.

        Returns:
            Integer: A random integer uniformly distributed in the range [min, max].

        """
        return self.thisptr.get().uniform_int(min, max)

    def gaussian(self, Real mean, Real sigma):
        """gaussian(mean, sigma) -> Real

        Return a Gaussian variate with the given mean and standard deviation.

        Args:
            mean (Real): The mean value.
            sigma (Real): The standard deviation.

        Returns:
            Real: A random number from a Gaussian distribution.

        """
        return self.thisptr.get().gaussian(mean, sigma)

    def binomial(self, Real p, Integer n):
        """binomial(p, n) -> Integer

        Return a random integer from the binomial distribution,
        the number of successes in n independent trials with probability p.

        Args:
            p (Real): A probability.
            n (Integer): The number of trials.

        Returns:
            Integer: A random integer from a binomial distribution.

        """
        return self.thisptr.get().binomial(p, n)

    def seed(self, val = None):
        """seed(val=None)

        Reset the random number seed.

        Args:
            val (Integer, optional): A new seed.
                If no seed is given, reset the seed by the current time.

        """
        if val is None:
            self.thisptr.get().seed()
        else:
            self.thisptr.get().seed(<Integer> val)

cdef GSLRandomNumberGenerator GSLRandomNumberGenerator_from_Cpp_RandomNumberGenerator(
    shared_ptr[Cpp_RandomNumberGenerator] rng):
    r = GSLRandomNumberGenerator()
    r.thisptr.swap(rng)
    return r
