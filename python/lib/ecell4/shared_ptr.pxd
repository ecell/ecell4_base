cdef extern from "<boost/shared_ptr.hpp>" namespace "boost":
    cdef cppclass shared_ptr[T]:
        shared_ptr()
        shared_ptr(T*)
        shared_ptr(shared_ptr[T]&)
        T* get()
        void swap(shared_ptr[T]&)

cdef extern from "<boost/pointer_cast.hpp>" namespace "boost":
    shared_ptr[T] dynamic_pointer_cast[T, U](shared_ptr[U]&)
    shared_ptr[T] static_pointer_cast[T, U](shared_ptr[U]&)
    shared_ptr[T] const_pointer_cast[T, U](shared_ptr[U]&)
