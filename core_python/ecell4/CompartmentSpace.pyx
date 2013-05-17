from cython.operator cimport dereference as deref


cdef class CompartmentSpaceVectorImpl:

    def __cinit__(self, Real volume):
        self.thisptr = new Cpp_CompartmentSpaceVectorImpl(volume)

    def __dealloc__(self):
        del self.thisptr

    def volume(self):
        return self.thisptr.volume()

    def num_species(self):
        return self.thisptr.num_species()

    def has_species(self, Species sp):
        return self.thisptr.has_species(deref(sp.thisptr))

    def num_molecules(self, Species sp):
        return self.thisptr.num_molecules(deref(sp.thisptr))

    def set_volume(self, Real volume):
        self.thisptr.set_volume(volume)

    def add_species(self, Species sp):
        self.thisptr.add_species(deref(sp.thisptr))

    def remove_species(self, Species sp):
        self.thisptr.remove_species(deref(sp.thisptr))

    def add_molecules(self, Species sp, Integer num):
        self.thisptr.add_molecules(deref(sp.thisptr), num)

    def remove_molecules(self, Species sp, Integer num):
        self.thisptr.remove_molecules(deref(sp.thisptr), num)
