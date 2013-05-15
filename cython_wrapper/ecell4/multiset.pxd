cdef extern from "<set>" namespace "std":
    cdef cppclass multiset[T]:
        multiset() except +
        multiset(multiset &) except+
        cppclass iterator:
            T& operator*()
            iterator operator++()
            iterator operator--()
            bint operator==(iterator)
            bint operator!=(iterator)
        iterator begin()
        iterator end()
