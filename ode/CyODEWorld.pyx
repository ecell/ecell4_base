# distutils: language = c++
# distutils: sources = ODEWorld.cpp

from libcpp.vector cimport vector

ctypedef double Real


cdef extern from "ODEWorld.hpp" namespace "ecell4::ode":
    cdef cppclass ODEWorld:
        ODEWorld(Real volume) except +
        Real & volume()
        void set_volume(Real & volume)
        int num_species()
#       bool has_species(Species sp)
#       void add_species(Species sp) 
#       void remove_species(Species sp)
#       Real num_molecules(Species sp)
#       void set_num_molecules(Species sp, Real const num)
#       void add_molecules(Species sp, Real num) 
#       void remove_molecules(Species sp, Real num)
        
        
cdef class PyOdeWorld:
    cdef ODEWorld *thisptr
    def __cinit__(self, Real vol):
        self.thisptr = new ODEWorld(vol)
    def __dealloc__(self):
        del self.thisptr
    def volume(self):
        return self.thisptr.volume()

    def set_volume(self, Real vol):
        self.thisptr.set_volume(vol)

    def num_species(self):
        return self.thisptr.num_species()

