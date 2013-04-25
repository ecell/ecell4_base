
from cython.operator cimport dereference as deref
from cython cimport address, declare
from libcpp.vector cimport vector
from libcpp.string cimport string
        

cdef class PyParticle:
    #cdef Position3 *thisptr
    def __cinit__(self, PySpecies sp, PyPosition3 pos, Real radius, Real D):
        self.thisptr = new Particle( 
                deref(sp.thisptr), 
                deref(pos.thisptr), 
                radius, D)
    def __dealloc__(self):
        del self.thisptr
    def position(self):
        cdef Position3 pos = self.thisptr.position()
        return to_PyPosition3( address(pos) )
    def radius(self):
        return self.this.radius()
    def D(self):
        return self.thisptr.D()
    def species(self):
        return to_PySpecies( address(self.thisptr.species()) )
