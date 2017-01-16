cdef class GSLRandomNumberGenerator:
    """A random number generator using the GNU Scientific Library (GSL).

    GSLRandomNumberGenerator(Integer myseed=None)

    """

    def __init__(self, val=None):
        """Constructor.

        Parameters
        ----------
        val : Integer or string, optional
            A seed for the random number generation (Integer),
            or a HDF5 filename (string).

        """
        pass

    def __cinit__(self, val=None):
        if val is None:
            self.thisptr = new shared_ptr[Cpp_RandomNumberGenerator](
                <Cpp_RandomNumberGenerator*> (new Cpp_GSLRandomNumberGenerator()))
        elif isinstance(val, str):
            self.thisptr = new shared_ptr[Cpp_RandomNumberGenerator](
                <Cpp_RandomNumberGenerator*> (
                    new Cpp_GSLRandomNumberGenerator(tostring(val))))
        else:
            self.thisptr = new shared_ptr[Cpp_RandomNumberGenerator](
                <Cpp_RandomNumberGenerator*> (
                    new Cpp_GSLRandomNumberGenerator(<Integer>val)))

    def __dealloc__(self):
        del self.thisptr

    def uniform(self, Real min, Real max):
        """uniform(min, max) -> Real

        Return a uniform random number within the given range.

        Parameters
        ----------
        min : Real
            The minimum value in the range.
        max : Real
            The maximum value in the range.

        Returns
        -------
        Real:
            A random number uniformly distributed in the range [min, max).

        """
        return self.thisptr.get().uniform(min, max)

    def uniform_int(self, Integer min, Integer max):
        """uniform_int(min, max) -> Integer

        Return a uniform random number within the given range.

        Parameters
        ----------
        min : Real
            The minimum value in the range.
        max : Real
            The maximum value in the range.

        Returns
        -------
        Integer:
            A random integer uniformly distributed in the range [min, max].

        """
        return self.thisptr.get().uniform_int(min, max)

    def gaussian(self, Real sigma, mean = None):
        """gaussian(sigma, mean = None) -> Real

        Return a Gaussian variate with the given mean and standard deviation.

        Parameters
        ----------
        sigma : Real
            The standard deviation.
        mean : Real
            The mean value.

        Returns
        -------
        Real:
            A random number from a Gaussian distribution.

        """
        if mean is None:
            return self.thisptr.get().gaussian(sigma)
        else:
            return self.thisptr.get().gaussian(sigma, <Real>mean)

    def binomial(self, Real p, Integer n):
        """binomial(p, n) -> Integer

        Return a random integer from the binomial distribution,
        the number of successes in n independent trials with probability p.

        Parameters
        ----------
        p : Real
            A probability.
        n : Integer
            The number of trials.

        Returns
        -------
        Integer:
            A random integer from a binomial distribution.

        """
        return self.thisptr.get().binomial(p, n)

    def seed(self, val = None):
        """seed(val=None)

        Reset the random number seed.

        Parameters
        ----------
        val : Integer, optional
            A new seed. If no seed is given, reset the seed by the current time.

        """
        if val is None:
            self.thisptr.get().seed()
        else:
            self.thisptr.get().seed(<Integer> val)

    def save(self, filename):
        """save(filename)

        Save the random number generator state to a file.

        Parameters
        ----------
        filename : str
            A filename to save to

        """
        self.thisptr.get().save(tostring(filename))

    def load(self, filename):
        """load(filename)

        Load the random number generator state from a file.

        Parameters
        ----------
        filename : str
            A filename to load from

        """
        self.thisptr.get().load(tostring(filename))

cdef GSLRandomNumberGenerator GSLRandomNumberGenerator_from_Cpp_RandomNumberGenerator(
    shared_ptr[Cpp_RandomNumberGenerator] rng):
    r = GSLRandomNumberGenerator()
    r.thisptr.swap(rng)
    return r
