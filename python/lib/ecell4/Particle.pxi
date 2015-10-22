from cython.operator cimport dereference as deref
from cython cimport address, declare
from libcpp.vector cimport vector
from libcpp.pair cimport pair
from libcpp.string cimport string
cimport util


cdef class ParticleID:
    """A class representing an ID of each particle

    ParticleID(value)

    """

    def __init__(self, value = None):
        """Constructor.

        Parameters
        ----------
        value : tuple
            A pair of integers, lot and serial.

        """
        pass

    def __cinit__(self, value = None):
        cdef pair[int, unsigned long long] val
        if value is None:
            self.thisptr = new Cpp_ParticleID()
        else:
            val.first = value[0]
            val.second = value[1]
            self.thisptr = new Cpp_ParticleID(val)

    def __dealloc__(self):
        del self.thisptr

    def __richcmp__(ParticleID self, ParticleID rhs, int op):
        cdef int compare
        if deref(self.thisptr) > deref(rhs.thisptr):
            compare = 1
        elif deref(self.thisptr) < deref(rhs.thisptr):
            compare = -1
        else: # self == rhs
            compare = 0
        return util.richcmp_helper(compare, op)

    def lot(self):
        """Return the first value."""
        return self.thisptr.lot()

    def serial(self):
        """Return the second value."""
        return self.thisptr.serial()

cdef ParticleID ParticleID_from_Cpp_ParticleID(Cpp_ParticleID* p):
    cdef Cpp_ParticleID *new_obj = new Cpp_ParticleID(<Cpp_ParticleID> deref(p))
    r = ParticleID((0, 0))
    del r.thisptr
    r.thisptr = new_obj
    return r

cdef class Particle:
    """A class representing a particle

    Particle(Species sp, Real3 pos, Real radius, Real D)

    """

    def __init__(self, Species sp, Real3 pos, Real radius, Real D):
        """Constructor.

        Parameters
        ----------
        sp : Species
            A species, which the particle belongs to.
        pos : Real3
            A position of the particle.
        radius : Real
            A radius of the particle. This must be positive.
        D : Real
            A diffusion constant. This must be positive.

        """
        pass

    def __cinit__(self, Species sp, Real3 pos, Real radius, Real D):
        self.thisptr = new Cpp_Particle(
            deref(sp.thisptr), deref(pos.thisptr), radius, D)

    def __dealloc__(self):
        del self.thisptr

    def position(self):
        """position() -> Real3

        Return the position.

        """
        cdef Cpp_Real3 pos = self.thisptr.position()
        return Real3_from_Cpp_Real3(address(pos))

    def radius(self):
        """Return the radius."""
        return self.thisptr.radius()

    def D(self):
        """Return the diffusion coefficient."""
        return self.thisptr.D()

    def species(self):
        """species() -> Species

        Return the species.

        """
        cdef Cpp_Species sp = self.thisptr.species()
        return Species_from_Cpp_Species(address(sp))

cdef Particle Particle_from_Cpp_Particle(Cpp_Particle* p):
    cdef Cpp_Particle *new_obj = new Cpp_Particle(<Cpp_Particle> deref(p))
    r = Particle(Species(), Real3(0, 0, 0), 0, 0)
    del r.thisptr
    r.thisptr = new_obj
    return r
