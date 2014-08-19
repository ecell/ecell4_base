cdef class Observer:

    def __cinit__(self):
        self.thisptr = new shared_ptr[Cpp_Observer](
            <Cpp_Observer*>(new Cpp_FixedIntervalNumberObserver(
                0.0, vector[string]()))) #XXX: DUMMY

    def __dealloc__(self):
        del self.thisptr

    def next_time(self):
        return self.thisptr.get().next_time()

cdef class FixedIntervalNumberObserver:

    def __cinit__(self, Real dt, vector[string] species):
        self.thisptr = new shared_ptr[Cpp_FixedIntervalNumberObserver](
            new Cpp_FixedIntervalNumberObserver(dt, species))

    def __dealloc__(self):
        del self.thisptr

    def next_time(self):
        return self.thisptr.get().next_time()

    def data(self):
        cdef vector[vector[Real]] d = self.thisptr.get().data()
        retval = []
        cdef vector[vector[Real]].iterator it = d.begin()
        while it != d.end():
            retval.append(deref(it))
            inc(it)
        return retval

    def as_base(self):
        retval = Observer()
        del retval.thisptr
        retval.thisptr = new shared_ptr[Cpp_Observer](
            <shared_ptr[Cpp_Observer]>deref(self.thisptr))
        return retval

cdef class NumberObserver:

    def __cinit__(self, vector[string] species):
        self.thisptr = new shared_ptr[Cpp_NumberObserver](
            new Cpp_NumberObserver(species))

    def __dealloc__(self):
        del self.thisptr

    def next_time(self):
        return self.thisptr.get().next_time()

    def data(self):
        cdef vector[vector[Real]] d = self.thisptr.get().data()
        retval = []
        cdef vector[vector[Real]].iterator it = d.begin()
        while it != d.end():
            retval.append(deref(it))
            inc(it)
        return retval

    def as_base(self):
        retval = Observer()
        del retval.thisptr
        retval.thisptr = new shared_ptr[Cpp_Observer](
            <shared_ptr[Cpp_Observer]>deref(self.thisptr))
        return retval
