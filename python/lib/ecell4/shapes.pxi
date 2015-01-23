from cython.operator cimport dereference as deref
from cython cimport address


cdef class Shape:

    def __cinit__(self):
        self.thisptr = new shared_ptr[Cpp_Shape](
            <Cpp_Shape*>(new Cpp_Sphere())) #XXX: DUMMY

    def __dealloc__(self):
        del self.thisptr

    def is_inside(self, Real3 pos):
        return self.thisptr.get().is_inside(deref(pos.thisptr))

    def dimension(self):
        return self.thisptr.get().dimension()

cdef class Sphere:

    def __cinit__(self, Real3 center, Real radius):
        self.thisptr = new shared_ptr[Cpp_Sphere](
            new Cpp_Sphere(deref(center.thisptr), radius))

    def __dealloc__(self):
        del self.thisptr

    def dimension(self):
        return self.thisptr.get().dimension()

    def distance(self, Real3 pos):
        return self.thisptr.get().distance(deref(pos.thisptr))

    def is_inside(self, Real3 pos):
        return self.thisptr.get().is_inside(deref(pos.thisptr))

    def surface(self):
        cdef Cpp_SphericalSurface shape = self.thisptr.get().surface()
        return SphericalSurface_from_Cpp_SphericalSurface(address(shape))

    def as_base(self):
        cdef shared_ptr[Cpp_Shape] *new_obj = new shared_ptr[Cpp_Shape](
            <Cpp_Shape*>(new Cpp_Sphere(<Cpp_Sphere> deref(self.thisptr.get()))))
        retval = Shape()
        del retval.thisptr
        retval.thisptr = new_obj
        return retval

cdef class SphericalSurface:

    def __cinit__(self, Real3 center, Real radius):
        self.thisptr = new shared_ptr[Cpp_SphericalSurface](
            new Cpp_SphericalSurface(deref(center.thisptr), radius))

    def __dealloc__(self):
        del self.thisptr

    def dimension(self):
        return self.thisptr.get().dimension()

    def distance(self, Real3 pos):
        return self.thisptr.get().distance(deref(pos.thisptr))

    def is_inside(self, Real3 pos):
        return self.thisptr.get().is_inside(deref(pos.thisptr))

    def inside(self):
        cdef Cpp_Sphere shape = self.thisptr.get().inside()
        return Sphere_from_Cpp_Sphere(address(shape))

    def as_base(self):
        cdef shared_ptr[Cpp_Shape] *new_obj = new shared_ptr[Cpp_Shape](
            <Cpp_Shape*>(new Cpp_SphericalSurface(<Cpp_SphericalSurface> deref(self.thisptr.get()))))
        retval = Shape()
        del retval.thisptr
        retval.thisptr = new_obj
        return retval

cdef class PlanarSurface:
    def __cinit__(self, Real3 origin, Real3 e0, Real3 e1):
        self.thisptr = new shared_ptr[Cpp_PlanarSurface](
            new Cpp_PlanarSurface(deref(origin.thisptr),
                deref(e0.thisptr), deref(e1.thisptr)))

    def __dealloc__(self):
        del self.thisptr

    def dimension(self):
        return self.thisptr.get().dimension()

    def is_inside(self, Real3 pos):
        return self.thisptr.get().is_inside(deref(pos.thisptr))

    def as_base(self):
        cdef shared_ptr[Cpp_Shape] *new_obj = new shared_ptr[Cpp_Shape](
            <Cpp_Shape*>(new Cpp_PlanarSurface(<Cpp_PlanarSurface> deref(self.thisptr.get()))))
        retval = Shape()
        del retval.thisptr
        retval.thisptr = new_obj
        return retval

cdef class AABB:

    def __cinit__(self, Real3 lower, Real3 upper):
        self.thisptr = new shared_ptr[Cpp_AABB](
            new Cpp_AABB(deref(lower.thisptr), deref(upper.thisptr)))

    def __dealloc__(self):
        del self.thisptr

    def dimension(self):
        return self.thisptr.get().dimension()

    def distance(self, Real3 pos):
        return self.thisptr.get().distance(deref(pos.thisptr))

    def is_inside(self, Real3 pos):
        return self.thisptr.get().is_inside(deref(pos.thisptr))

    def upper(self):
        cdef Cpp_Real3 pos = self.thisptr.get().upper()
        return Real3_from_Cpp_Real3(address(pos))

    def lower(self):
        cdef Cpp_Real3 pos = self.thisptr.get().lower()
        return Real3_from_Cpp_Real3(address(pos))

    def as_base(self):
        cdef shared_ptr[Cpp_Shape] *new_obj = new shared_ptr[Cpp_Shape](
            <Cpp_Shape*>(new Cpp_AABB(<Cpp_AABB> deref(self.thisptr.get()))))
        retval = Shape()
        del retval.thisptr
        retval.thisptr = new_obj
        return retval

cdef Sphere Sphere_from_Cpp_Sphere(Cpp_Sphere* shape):
    cdef shared_ptr[Cpp_Sphere] *new_obj = new shared_ptr[Cpp_Sphere](
        new Cpp_Sphere(<Cpp_Sphere> deref(shape)))
    retval = Sphere(Real3(0, 0, 0), 0)
    del retval.thisptr
    retval.thisptr = new_obj
    return retval

cdef SphericalSurface SphericalSurface_from_Cpp_SphericalSurface(Cpp_SphericalSurface* shape):
    cdef shared_ptr[Cpp_SphericalSurface] *new_obj = new shared_ptr[Cpp_SphericalSurface](
        new Cpp_SphericalSurface(<Cpp_SphericalSurface> deref(shape)))
    retval = SphericalSurface(Real3(0, 0, 0), 0)
    del retval.thisptr
    retval.thisptr = new_obj
    return retval

cdef AABB AABB_from_Cpp_AABB(Cpp_AABB* shape):
    cdef shared_ptr[Cpp_AABB] *new_obj = new shared_ptr[Cpp_AABB](
        new Cpp_AABB(<Cpp_AABB> deref(shape)))
    retval = AABB(Real3(0, 0, 0), Real3(0, 0, 0))
    del retval.thisptr
    retval.thisptr = new_obj
    return retval

