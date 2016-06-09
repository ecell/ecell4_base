from cython.operator cimport dereference as deref, preincrement as inc
from cython cimport address

from libcpp.pair cimport pair
from libcpp.vector cimport vector
from libcpp.string cimport string
cimport util

cdef class PlanarSurfaceID:
    """A class representing an ID of each PlanarSurface(PlanarSurface):
    PlanarSurfaceID(value)
    """
    def __cinit__(self, value=None):
        cdef pair[int, unsigned long long] val
        if value is None:
            self.thisptr = new Cpp_PlanarSurfaceID()
        else:
            val.first = value[0]
            val.second = value[1]
            self.thisptr = new Cpp_PlanarSurfaceID(val)
    def __dealloc__(self):
        del self.thisptr

    def __richcmp__(PlanarSurfaceID self, PlanarSurfaceID rhs, int op):
        cdef int compare
        if deref(self.thisptr) > deref(rhs.thisptr):
            compare = 1
        elif deref(self.thisptr) < deref(rhs.thisptr):
            compare = -1
        else:
            compare = 0
        return util.richcmp_helper(compare, op)

    def lot(self):
        """Return the first value."""
        return self.thisptr.lot()

    def serial(self):
        """Return the second value."""
        return self.thisptr.serial()

    def __reduce__(self):
        return (PlanarSurfaceID, ((self.lot(), self.serial()), ))

cdef PlanarSurfaceID PlanarSurfaceID_from_Cpp_PlanarSurfaceID(Cpp_PlanarSurfaceID *p):
    cdef Cpp_PlanarSurfaceID *new_obj = new Cpp_PlanarSurfaceID(<Cpp_PlanarSurfaceID> deref(p))
    r = PlanarSurfaceID((0, 0))
    del r.thisptr
    r.thisptr = new_obj
    return r





