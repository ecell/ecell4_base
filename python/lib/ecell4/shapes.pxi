from cython.operator cimport dereference as deref
from cython cimport address


cdef class Shape:

    def __cinit__(self):
        self.thisptr = <Cpp_Shape*>(
            new Cpp_Sphere(Cpp_Position3(0.0, 0.0, 0.0), 0.0)) #XXX: DUMMY

    def __dealloc__(self):
        del self.thisptr

    def is_inside(self, Position3 pos):
        return self.thisptr.is_inside(deref(pos.thisptr))

    def dimension(self):
        return self.thisptr.dimension()

cdef class Sphere:

    def __cinit__(self, Position3 center, Real radius):
        self.thisptr = new Cpp_Sphere(deref(center.thisptr), radius)

    def __dealloc__(self):
        del self.thisptr

    def dimension(self):
        return self.thisptr.dimension()

    def distance(self, Position3 pos):
        return self.thisptr.distance(deref(pos.thisptr))

    def is_inside(self, Position3 pos):
        return self.thisptr.is_inside(deref(pos.thisptr))

    def surface(self):
        cdef Cpp_SphericalSurface shape = self.thisptr.surface()
        return SphericalSurface_from_Cpp_SphericalSurface(address(shape))

    def as_base(self):
        cdef Cpp_Shape *new_obj = <Cpp_Shape*>(new Cpp_Sphere(<Cpp_Sphere> deref(self.thisptr)))
        retval = Shape()
        del retval.thisptr
        retval.thisptr = new_obj
        return retval

cdef class SphericalSurface:

    def __cinit__(self, Position3 center, Real radius):
        self.thisptr = new Cpp_SphericalSurface(deref(center.thisptr), radius)

    def __dealloc__(self):
        del self.thisptr

    def dimension(self):
        return self.thisptr.dimension()

    def distance(self, Position3 pos):
        return self.thisptr.distance(deref(pos.thisptr))

    def is_inside(self, Position3 pos):
        return self.thisptr.is_inside(deref(pos.thisptr))

    def inside(self):
        cdef Cpp_Sphere shape = self.thisptr.inside()
        return Sphere_from_Cpp_Sphere(address(shape))

    def as_base(self):
        cdef Cpp_Shape *new_obj = <Cpp_Shape*>(
            new Cpp_SphericalSurface(<Cpp_SphericalSurface> deref(self.thisptr)))
        retval = Shape()
        del retval.thisptr
        retval.thisptr = new_obj
        return retval

cdef Sphere Sphere_from_Cpp_Sphere(Cpp_Sphere* shape):
    cdef Cpp_Sphere *new_obj = new Cpp_Sphere(<Cpp_Sphere> deref(shape))
    retval = Sphere(Position3(0, 0, 0), 0)
    del retval.thisptr
    retval.thisptr = new_obj
    return retval

cdef SphericalSurface SphericalSurface_from_Cpp_SphericalSurface(Cpp_SphericalSurface* shape):
    cdef Cpp_SphericalSurface *new_obj = new Cpp_SphericalSurface(<Cpp_SphericalSurface> deref(shape))
    retval = SphericalSurface(Position3(0, 0, 0), 0)
    del retval.thisptr
    retval.thisptr = new_obj
    return retval
