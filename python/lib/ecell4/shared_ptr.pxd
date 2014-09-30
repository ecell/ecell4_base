cdef extern from "<boost/shared_ptr.hpp>" namespace "boost":
    cdef cppclass shared_ptr[T]:
        shared_ptr(T *ptr)
        shared_ptr(shared_ptr[T] ptr)
        T* get()
        void swap(shared_ptr[T]&)
