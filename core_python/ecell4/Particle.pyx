from cython.operator cimport dereference as deref
from cython cimport address, declare
from libcpp.vector cimport vector
from libcpp.pair cimport pair
from libcpp.string cimport string


cdef class ParticleID:

    def __cinit__(self, lot, serial):
        cdef pair[int, unsigned long long] val
        val.first = lot
        val.second = serial
        self.thisptr = new Cpp_ParticleID(val)

    def __dealloc__(self):
        del self.thisptr

    def lot(self):
        return self.thisptr.lot()

    def serial(self):
        return self.thisptr.serial()

cdef class Particle:

    def __cinit__(self, Species sp, Position3 pos, Real radius, Real D):
        self.thisptr = new Cpp_Particle(
            deref(sp.thisptr), deref(pos.thisptr), radius, D)

    def __dealloc__(self):
        del self.thisptr

    def position(self):
        cdef Cpp_Position3 pos = self.thisptr.position()
        return Position3_from_Cpp_Position3(address(pos))

    def radius(self):
        return self.thisptr.radius()

    def D(self):
        return self.thisptr.D()

    def species(self):
        return Species_from_Cpp_Species(address(self.thisptr.species()))

cdef ParticleID ParticleID_from_Cpp_ParticleID(Cpp_ParticleID* p):
    cdef Cpp_ParticleID *new_obj = new Cpp_ParticleID(<Cpp_ParticleID> deref(p))
    r = ParticleID(0, 0)
    del r.thisptr
    r.thisptr = new_obj
    return r

cdef Particle Particle_from_Cpp_Particle(Cpp_Particle* p):
    cdef Cpp_Particle *new_obj = new Cpp_Particle(<Cpp_Particle> deref(p))
    r = Particle(Species(), Position3(0, 0, 0), 0, 0)
    del r.thisptr
    r.thisptr = new_obj
    return r
