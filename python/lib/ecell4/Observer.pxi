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

    def num_steps(self):
        return self.thisptr.get().num_steps()

    def data(self):
        cdef vector[vector[Real]] d = self.thisptr.get().data()
        retval = []
        cdef vector[vector[Real]].iterator it = d.begin()
        while it != d.end():
            retval.append(deref(it))
            inc(it)
        return retval

    def targets(self):
        cdef vector[Cpp_Species] species = self.thisptr.get().targets()

        retval = []
        cdef vector[Cpp_Species].iterator it = species.begin()
        while it != species.end():
            retval.append(
                 Species_from_Cpp_Species(
                     <Cpp_Species*>(address(deref(it)))))
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

    def targets(self):
        cdef vector[Cpp_Species] species = self.thisptr.get().targets()

        retval = []
        cdef vector[Cpp_Species].iterator it = species.begin()
        while it != species.end():
            retval.append(
                 Species_from_Cpp_Species(
                     <Cpp_Species*>(address(deref(it)))))
            inc(it)
        return retval

    def as_base(self):
        retval = Observer()
        del retval.thisptr
        retval.thisptr = new shared_ptr[Cpp_Observer](
            <shared_ptr[Cpp_Observer]>deref(self.thisptr))
        return retval

cdef class FixedIntervalHDF5Observer:

    def __cinit__(self, Real dt, string filename):
        self.thisptr = new shared_ptr[Cpp_FixedIntervalHDF5Observer](
            new Cpp_FixedIntervalHDF5Observer(dt, filename))

    def __dealloc__(self):
        del self.thisptr

    def next_time(self):
        return self.thisptr.get().next_time()

    def num_steps(self):
        return self.thisptr.get().num_steps()

    def filename(self):
        return self.thisptr.get().filename()

    def as_base(self):
        retval = Observer()
        del retval.thisptr
        retval.thisptr = new shared_ptr[Cpp_Observer](
            <shared_ptr[Cpp_Observer]>deref(self.thisptr))
        return retval
