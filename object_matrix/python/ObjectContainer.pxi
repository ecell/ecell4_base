include "Sphere.pxi"

cdef extern from "_object_container.hpp":
    ctypedef struct __impl_object_container "my_object_container_type":
        void erase(int key)
        int size()
        int matrix_size()
        double cell_size()


    ctypedef struct __impl_sphere_ref "sphere_ref":
        double x()
        double y()
        double z()
        double radius()

    __impl_sphere *__impl_sphere_ref_clone(__impl_sphere_ref *)

    void __impl_sphere_ref_del "delete" (__impl_sphere_ref *)

    ctypedef struct __impl_object_container_iterneighbors_gen_value "take_neighbor_collector::value_type":
        __impl_sphere_ref *first
        double second

    ctypedef struct __impl_object_container_iterneighbors_gen "take_neighbor_collector::generator_type":
        __impl_object_container_iterneighbors_gen_value* current "operator*"()
        __impl_object_container_iterneighbors_gen_value* next "operator++"()

    __impl_object_container_iterneighbors_gen *__impl_object_container_iterneighbors(__impl_object_container *, __impl_sphere *)

    __impl_object_container_iterneighbors_gen *__impl_object_container_iterneighbors_cyclic(__impl_object_container *, __impl_sphere *)

    __impl_object_container*__impl_object_container_new "new my_object_container_type" (double world_size, int cells_per_side)

    void __impl_object_container_del "delete" (__impl_object_container *)

    void __impl_object_container_iterneighbors_gen_del "delete" (__impl_object_container_iterneighbors_gen *)

    int __impl_object_container_insert(__impl_object_container *,
            int key, __impl_sphere *)

    __impl_sphere_ref *__impl_object_container_find(__impl_object_container *,
            int key)
