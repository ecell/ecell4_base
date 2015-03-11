cdef class Observer:

    def __cinit__(self):
        self.thisptr = new shared_ptr[Cpp_Observer](
            <Cpp_Observer*>(new Cpp_FixedIntervalNumberObserver(
                0.0, vector[string]()))) #XXX: DUMMY

    def __dealloc__(self):
        del self.thisptr

    def next_time(self):
        return self.thisptr.get().next_time()

    def reset(self):
        self.thisptr.get().reset()

cdef class FixedIntervalNumberObserver:

    def __cinit__(self, Real dt, species):
        cdef vector[string] cpp_species
        for serial in species:
            cpp_species.push_back(tostring(serial))
        self.thisptr = new shared_ptr[Cpp_FixedIntervalNumberObserver](
            new Cpp_FixedIntervalNumberObserver(dt, cpp_species))

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

    def reset(self):
        self.thisptr.get().reset()

cdef class NumberObserver:

    def __cinit__(self, species):
        cdef vector[string] cpp_species
        for serial in species:
            cpp_species.push_back(tostring(serial))
        self.thisptr = new shared_ptr[Cpp_NumberObserver](
            new Cpp_NumberObserver(cpp_species))

    def __dealloc__(self):
        del self.thisptr

    def num_steps(self):
        return self.thisptr.get().num_steps()

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

    def reset(self):
        self.thisptr.get().reset()

cdef class FixedIntervalHDF5Observer:

    def __cinit__(self, Real dt, filename):
        self.thisptr = new shared_ptr[Cpp_FixedIntervalHDF5Observer](
            new Cpp_FixedIntervalHDF5Observer(dt, tostring(filename)))

    def __dealloc__(self):
        del self.thisptr

    def next_time(self):
        return self.thisptr.get().next_time()

    def num_steps(self):
        return self.thisptr.get().num_steps()

    def filename(self):
        return self.thisptr.get().filename().decode('UTF-8')

    def as_base(self):
        retval = Observer()
        del retval.thisptr
        retval.thisptr = new shared_ptr[Cpp_Observer](
            <shared_ptr[Cpp_Observer]>deref(self.thisptr))
        return retval

    def reset(self):
        self.thisptr.get().reset()

cdef class FixedIntervalCSVObserver:

    def __cinit__(self, Real dt, filename):
        self.thisptr = new shared_ptr[Cpp_FixedIntervalCSVObserver](
            new Cpp_FixedIntervalCSVObserver(dt, tostring(filename)))

    def __dealloc__(self):
        del self.thisptr

    def next_time(self):
        return self.thisptr.get().next_time()

    def num_steps(self):
        return self.thisptr.get().num_steps()

    def log(self, w):
        cdef Space space = w.as_base()
        self.thisptr.get().log(space.thisptr.get())

    def filename(self):
        return self.thisptr.get().filename().decode('UTF-8')

    def as_base(self):
        retval = Observer()
        del retval.thisptr
        retval.thisptr = new shared_ptr[Cpp_Observer](
            <shared_ptr[Cpp_Observer]>deref(self.thisptr))
        return retval

    def reset(self):
        self.thisptr.get().reset()

cdef class FixedIntervalTrajectoryObserver:

    def __cinit__(self, Real dt, pids, resolve_boundary=None):
        cdef vector[Cpp_ParticleID] tmp
        for pid in pids:
            tmp.push_back(deref((<ParticleID>pid).thisptr))
        if resolve_boundary is None:
            self.thisptr = new shared_ptr[Cpp_FixedIntervalTrajectoryObserver](
                new Cpp_FixedIntervalTrajectoryObserver(dt, tmp))
        else:
            self.thisptr = new shared_ptr[Cpp_FixedIntervalTrajectoryObserver](
                new Cpp_FixedIntervalTrajectoryObserver(dt, tmp, <bool>resolve_boundary))

    def __dealloc__(self):
        del self.thisptr

    def next_time(self):
        return self.thisptr.get().next_time()

    def num_steps(self):
        return self.thisptr.get().num_steps()

    def data(self):
        cdef vector[vector[Cpp_Real3]] d = self.thisptr.get().data()
        retval = []
        cdef vector[vector[Cpp_Real3]].iterator it = d.begin()
        cdef vector[Cpp_Real3].iterator it2
        while it != d.end():
            it2 = deref(it).begin()
            retval.append([])
            while it2 != deref(it).end():
                retval[-1].append(Real3_from_Cpp_Real3(address(deref(it2))))
                inc(it2)
            inc(it)
        return retval

    def as_base(self):
        retval = Observer()
        del retval.thisptr
        retval.thisptr = new shared_ptr[Cpp_Observer](
            <shared_ptr[Cpp_Observer]>deref(self.thisptr))
        return retval

    def reset(self):
        self.thisptr.get().reset()

cdef class TimingNumberObserver:

    def __cinit__(self, vector[double] t, species):  #XXX: vector[Real]
        cdef vector[string] cpp_species
        for serial in species:
            cpp_species.push_back(tostring(serial))
        self.thisptr = new shared_ptr[Cpp_TimingNumberObserver](
            new Cpp_TimingNumberObserver(t, cpp_species))

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

    def reset(self):
        self.thisptr.get().reset()
