from cython.operator cimport dereference as deref


cdef class Shape:

    def __cinit__(self):
        self.thisptr = new shared_ptr[Cpp_Shape](
            <Cpp_Shape*>(new Cpp_Sphere(
                Cpp_Position3(0.0, 0.0, 0.0), 0.0))) #XXX: DUMMY

    def __dealloc__(self):
        del self.thisptr

    def is_inside(self, Position3 pos):
        return self.thisptr.get().is_inside(deref(pos.thisptr))

cdef class Sphere:

    def __cinit__(self, Position3 center, Real radius):
        self.thisptr = new shared_ptr[Cpp_Sphere](
            new Cpp_Sphere(deref(center.thisptr), radius))

    def __dealloc__(self):
        del self.thisptr

    def distance(self, Position3 pos):
        return self.thisptr.get().distance(deref(pos.thisptr))

    def is_inside(self, Position3 pos):
        return self.thisptr.get().is_inside(deref(pos.thisptr))

    def as_base(self):
        retval = Shape()
        del retval.thisptr
        retval.thisptr = new shared_ptr[Cpp_Shape](
            <shared_ptr[Cpp_Shape]>deref(self.thisptr))
        return retval
