from cython.operator cimport dereference as deref
from cython cimport address, declare
from libcpp.vector cimport vector
from libcpp.string cimport string


cdef class Particle:

    def __cinit__(self, Species sp, Position3 pos, Real radius, Real D):
        self.thisptr = new Cpp_Particle(
            deref(sp.thisptr), deref(pos.thisptr), radius, D)

    def __dealloc__(self):
        del self.thisptr

    def position(self):
        cdef Cpp_Position3 pos = self.thisptr.position()
        return Cpp_Position3_to_Position3(address(pos))

    def radius(self):
        return self.this.radius()

    def D(self):
        return self.thisptr.D()

    def species(self):
        return Cpp_Species_to_Species(address(self.thisptr.species()))
