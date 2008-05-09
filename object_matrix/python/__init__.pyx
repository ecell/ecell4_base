include "python_api.pxi"
include "pyrex_helper.pxi"
include "Position.pxi"
include "Sphere.pxi"
include "ObjectContainer.pxi"

cdef class Position:
    cdef __impl_position *pimpl

    def __cinit__(self, double x = 0, double y = 0, double z = 0):
        self.pimpl = __impl_position_new(x, y, z)

    def __dealloc__(self):
        __impl_position_del(self.pimpl)

    property x:
        def __get__(self):
            return self.pimpl.x()
        def __set__(self, val):
            __helper_set_double(self.pimpl.x(), val)

    property y:
        def __get__(self):
            return self.pimpl.y()
        def __set__(self, val):
            __helper_set_double(self.pimpl.y(), val)

    property z:
        def __get__(self):
            return self.pimpl.z()
        def __set__(self, val):
            __helper_set_double(self.pimpl.z(), val)

cdef class Sphere:
    cdef __impl_sphere *pimpl

    def __cinit__(self, double x = 0, double y = 0, double z = 0, double radius = 0):
        self.pimpl = __impl_sphere_new(x, y, z, radius)

    def __dealloc__(self):
        __impl_sphere_del(self.pimpl)

    property x:
        def __get__(self):
            return self.pimpl.x()
        def __set__(self, val):
            __helper_set_double(self.pimpl.x(), val)

    property y:
        def __get__(self):
            return self.pimpl.y()
        def __set__(self, val):
            __helper_set_double(self.pimpl.y(), val)

    property z:
        def __get__(self):
            return self.pimpl.z()
        def __set__(self, val):
            __helper_set_double(self.pimpl.z(), val)

    property radius:
        def __get__(self):
            return self.pimpl.radius
        def __set__(self, val):
            __helper_set_double(self.pimpl.radius, val)

cdef class SphereRef:
    cdef __impl_sphere_ref *pimpl

    def __cinit__(self, int pimpl):
        self.pimpl = <__impl_sphere_ref*>pimpl

    def __repr__(self):
        return __helper_pystr_from_repr(self.pimpl)

    def __dealloc__(self):
        __impl_sphere_ref_del(self.pimpl)

    property x:
        def __get__(self):
            return self.pimpl.x()
        def __set__(self, val):
            __helper_set_double(self.pimpl.x(), val)

    property y:
        def __get__(self):
            return self.pimpl.y()
        def __set__(self, val):
            __helper_set_double(self.pimpl.y(), val)

    property z:
        def __get__(self):
            return self.pimpl.z()
        def __set__(self, val):
            __helper_set_double(self.pimpl.z(), val)

    property radius:
        def __get__(self):
            return self.pimpl.radius()
        def __set__(self, val):
            __helper_set_double(self.pimpl.radius(), val)

cdef class ObjectContainer

cdef class __ObjectContainer_iterneighbors_iter:
    cdef __impl_object_container_iterneighbors_gen *pimpl
    cdef __impl_object_container *oc_pimpl
    cdef __impl_sphere *sphere_pimpl

    def __cinit__(self, ObjectContainer oc, sphere):
        cdef __impl_sphere *sphere_pimpl
        if isinstance(sphere, Sphere):
            sphere_pimpl = __impl_sphere_clone((<Sphere>sphere).pimpl)
        elif isinstance(sphere, SphereRef):
            sphere_pimpl = __impl_sphere_ref_clone((<SphereRef>sphere).pimpl)
        else:
            raise ValueError("unsupported argument type")
        self.sphere_pimpl = sphere_pimpl
        self.oc_pimpl = oc.pimpl
        self.pimpl = NULL

    def __dealloc__(self):
        __impl_object_container_iterneighbors_gen_del(self.pimpl)
        __impl_sphere_del(self.sphere_pimpl)

    def next(self):
        cdef __impl_object_container_iterneighbors_gen_value *value
        if self.pimpl == NULL:
            self.pimpl = __impl_object_container_iterneighbors(
                    self.oc_pimpl, self.sphere_pimpl)
        value = <__impl_object_container_iterneighbors_gen_value *>self.pimpl.current()
        if value == NULL:
            return None
        retval = (SphereRef(<int>value.first), value.second)
        self.pimpl.next()
        return retval

cdef class __ObjectContainer_iterneighbors_cyclic_iter:
    cdef __impl_object_container_iterneighbors_gen *pimpl
    cdef __impl_object_container *oc_pimpl
    cdef __impl_sphere *sphere_pimpl

    def __cinit__(self, ObjectContainer oc, sphere):
        cdef __impl_sphere *sphere_pimpl
        if isinstance(sphere, Sphere):
            sphere_pimpl = __impl_sphere_clone((<Sphere>sphere).pimpl)
        elif isinstance(sphere, SphereRef):
            sphere_pimpl = __impl_sphere_ref_clone((<SphereRef>sphere).pimpl)
        else:
            raise ValueError("unsupported argument type")
        self.sphere_pimpl = sphere_pimpl
        self.oc_pimpl = oc.pimpl
        self.pimpl = NULL

    def __dealloc__(self):
        __impl_object_container_iterneighbors_gen_del(self.pimpl)
        __impl_sphere_del(self.sphere_pimpl)

    def next(self):
        cdef __impl_object_container_iterneighbors_gen_value *value
        if self.pimpl == NULL:
            self.pimpl = __impl_object_container_iterneighbors_cyclic(
                    self.oc_pimpl, self.sphere_pimpl)
        value = <__impl_object_container_iterneighbors_gen_value *>self.pimpl.current()
        if value == NULL:
            return None
        retval = (SphereRef(<int>value.first), value.second)
        self.pimpl.next()
        return retval

cdef class ObjectContainer:
    cdef __impl_object_container *pimpl

    def __cinit__(self, world_size = 1.0, cells_per_side = 1):
        self.pimpl = __impl_object_container_new(world_size, cells_per_side)

    def __dealloc__(self):
        __impl_object_container_del(self.pimpl)

    def __getitem__(self, key):
        s = <int>__impl_object_container_find(self.pimpl, key)
        if s == 0:
            return None
        return SphereRef(s)

    def __setitem__(self, key, Sphere sphere):
        __impl_object_container_insert(self.pimpl, key, sphere.pimpl)

    def __len__(self):
        return self.pimpl.size()


    def iterneighbors(self, sphere):
        return PyCallIter_New(__ObjectContainer_iterneighbors_iter(
                self, sphere).next, None)

    def iterneighbors_cyclic(self, sphere):
        return PyCallIter_New(__ObjectContainer_iterneighbors_cyclic_iter(
                self, sphere).next, None)

    property matrix_size:
        def __get__(self): return self.pimpl.matrix_size()

    property cell_size:
        def __get__(self): return self.pimpl.cell_size()

