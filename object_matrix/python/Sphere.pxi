cdef extern from "_sphere.hpp":
    ctypedef struct __impl_sphere "sphere<double>":
        double x "position.x"()
        double y "position.y"()
        double z "position.z"()
        double radius

    __impl_sphere *__impl_sphere_new(double x, double y, double z, double radius)

    __impl_sphere *__impl_sphere_clone(__impl_sphere *)

    void __impl_sphere_del "delete" (__impl_sphere *)
