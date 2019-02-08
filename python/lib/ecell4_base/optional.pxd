from libcpp cimport bool

cdef extern from "<boost/optional.hpp>" namespace "boost":
    cdef cppclass optional[T]:
        optional()
        optional(optional[T]&)
        T& get()
        void swap(optional[T]&)
        bool is_initialized()
