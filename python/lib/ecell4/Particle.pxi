from cython.operator cimport dereference as deref
from cython cimport address, declare
from libcpp.vector cimport vector
from libcpp.pair cimport pair
from libcpp.string cimport string
cimport util


cdef class ParticleID:

    def __cinit__(self, value = None):
        cdef pair[int, unsigned long long] val
        if value is None:
            self.thisptr = new Cpp_ParticleID()
        else:
            val.first = value[0]
            val.second = value[1]
            self.thisptr = new Cpp_ParticleID(val)

    # def __cinit__(self, lot, serial):
    #     cdef pair[int, unsigned long long] val
    #     val.first = lot
    #     val.second = serial
    #     self.thisptr = new Cpp_ParticleID(val)

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
        return self.thisptr.lot()

    def serial(self):
        return self.thisptr.serial()

cdef class Particle:

    def __cinit__(self, Species sp, Real3 pos, Real radius, Real D):
        self.thisptr = new Cpp_Particle(
            deref(sp.thisptr), deref(pos.thisptr), radius, D)

    def __dealloc__(self):
        del self.thisptr

    def position(self):
        cdef Cpp_Real3 pos = self.thisptr.position()
        return Real3_from_Cpp_Real3(address(pos))

    def radius(self):
        return self.thisptr.radius()

    def D(self):
        return self.thisptr.D()

    def species(self):
        cdef Cpp_Species sp = self.thisptr.species()
        return Species_from_Cpp_Species(address(sp))

cdef ParticleID ParticleID_from_Cpp_ParticleID(Cpp_ParticleID* p):
    cdef Cpp_ParticleID *new_obj = new Cpp_ParticleID(<Cpp_ParticleID> deref(p))
    r = ParticleID((0, 0))
    del r.thisptr
    r.thisptr = new_obj
    return r

cdef Particle Particle_from_Cpp_Particle(Cpp_Particle* p):
    cdef Cpp_Particle *new_obj = new Cpp_Particle(<Cpp_Particle> deref(p))
    r = Particle(Species(), Real3(0, 0, 0), 0, 0)
    del r.thisptr
    r.thisptr = new_obj
    return r
