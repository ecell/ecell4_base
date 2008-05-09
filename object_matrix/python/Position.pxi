cdef extern from "position.hpp":
    ctypedef struct __impl_position "position<double>":
        double x()
        double y()
        double z()

    __impl_position *__impl_position_new "new position<double>" (double x, double y, double z)

    void __impl_position_del "delete" (__impl_position *)
