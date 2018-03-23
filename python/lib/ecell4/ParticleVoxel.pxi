from cython.operator cimport dereference as deref
from cython cimport address, declare
from libcpp.vector cimport vector
from libcpp.pair cimport pair
from libcpp.string cimport string


cdef class ParticleVoxel:
    """A class representing a voxel in LatticeSpace.

    ParticleVoxel(Species sp, Integer coord, Real radius, Real D, loc=None)

    """

    def __init__(self, Species sp, Integer coord, Real radius, Real D, loc=None):
        """Constructor.

        Parameters
        ----------
        sp : Species
            The species.
        coord : Integer
            The coordinate given as an Integer.
        radius : Real
            The radius of a molecule.
        D : Real
            The diffusion rate of a molecule.
        loc : str, optional
            The location of a molecule.

        """
        pass  # XXX: Only used for doc string

    def __cinit__(self, Species sp, Integer coord, Real radius, Real D, loc=None):
        if loc is None:
            self.thisptr = new Cpp_ParticleVoxel(
                deref(sp.thisptr), coord, radius, D)
        else:
            self.thisptr = new Cpp_ParticleVoxel(
                deref(sp.thisptr), coord, radius, D, tostring(loc))

    def __dealloc__(self):
        del self.thisptr

    def coordinate(self):
        """Return the coordinate."""
        return self.thisptr.coordinate

    def D(self):
        """Return the diffusion coefficient."""
        return self.thisptr.D

    def radius(self):
        """Return the radius."""
        return self.thisptr.radius

    def species(self):
        """species() -> Species

        Return the species.

        """
        return Species_from_Cpp_Species(address(self.thisptr.species))

    def loc(self):
        """loc() -> str

        Return the location information as a string.

        """
        return self.thisptr.loc.decode('UTF-8')

    def __reduce__(self):
        return (ParticleVoxel, (self.species(), self.coordinate(), self.radius(), self.D(), self.loc()))

cdef ParticleVoxel Voxel_from_Cpp_Voxel(Cpp_ParticleVoxel* p):
    cdef Cpp_ParticleVoxel *new_obj = new Cpp_ParticleVoxel(<Cpp_ParticleVoxel> deref(p))
    r = ParticleVoxel(Species(), 0, 0, 0)
    del r.thisptr
    r.thisptr = new_obj
    return r
