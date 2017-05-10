import types

cdef class Observer:
    """A wrapper for a base class of Observers.

    Warning: This is mainly for developers.
    Do not use this for your simulation.

    """

    def __init__(self):
        """Constructor."""
        pass

    def __cinit__(self):
        self.thisptr = new shared_ptr[Cpp_Observer](
            <Cpp_Observer*>(new Cpp_FixedIntervalNumberObserver(
                0.0, vector[string]()))) #XXX: DUMMY

    def __dealloc__(self):
        del self.thisptr

    def next_time(self):
        """Return the next time for logging."""
        return self.thisptr.get().next_time()

    def reset(self):
        """Reset the internal state."""
        self.thisptr.get().reset()

cdef class FixedIntervalNumberObserver:
    """An ``Observer``class to log the number of molecules with the fixed
    step interval.
    This ``Observer`` logs at the current time first, and then keeps logging
    every after the interval.

    FixedIntervalNumberObserver(dt, species)

    """

    def __init__(self, Real dt, species):
        """Constructor.

        Parameters
        ----------
        dt : float
            A step interval for logging.
        species : list
            A list of strings, but not of ``Species``.
            The strings suggest serials of ``Species`` to be observed.

        """
        pass  # XXX: Only used for doc string

    def __cinit__(self, Real dt, species):
        cdef vector[string] cpp_species
        for serial in species:
            cpp_species.push_back(tostring(serial))
        self.thisptr = new shared_ptr[Cpp_FixedIntervalNumberObserver](
            new Cpp_FixedIntervalNumberObserver(dt, cpp_species))

    def __dealloc__(self):
        del self.thisptr

    def next_time(self):
        """Return the next time for logging."""
        return self.thisptr.get().next_time()

    def num_steps(self):
        """Return the number of steps."""
        return self.thisptr.get().num_steps()

    def data(self):
        """Return a list of the number of molecules you specified.

        Returns
        -------
        list:
            A list of lists of the numbers of molecules.
            The size of a return value is equal to ``num_steps``.
            Each element of a return value is a list consisting of
            time and the number of molecules specified at the construction.

        """
        cdef vector[vector[Real]] d = self.thisptr.get().data()
        retval = []
        cdef vector[vector[Real]].iterator it = d.begin()
        while it != d.end():
            retval.append(deref(it))
            inc(it)
        return retval

    def targets(self):
        """Return a list of ``Species``, which this ``Observer`` observes

        Returns
        -------
        list:
            A list of ``Species``. This is generated from arguments
            you gave at the construction.

        """
        cdef vector[Cpp_Species] species = self.thisptr.get().targets()

        retval = []
        cdef vector[Cpp_Species].iterator it = species.begin()
        while it != species.end():
            retval.append(
                 Species_from_Cpp_Species(
                     <Cpp_Species*>(address(deref(it)))))
            inc(it)
        return retval

    def save(self, filename):
        """Save data to an output with the given filename."""
        self.thisptr.get().save(tostring(filename))

    def as_base(self):
        """Clone self as a base class. This function is for developers."""
        retval = Observer()
        del retval.thisptr
        retval.thisptr = new shared_ptr[Cpp_Observer](
            <shared_ptr[Cpp_Observer]>deref(self.thisptr))
        return retval

    def reset(self):
        """Reset the internal state."""
        self.thisptr.get().reset()

cdef class NumberObserver:
    """An ``Observer``class to log the number of molecules.
    This ``Observer`` logs at the current time first, and then keeps logging
    every after simulation steps.
    Warning: This doesn't work with ODESimulator.

    NumberObserver(species)

    """

    def __init__(self, species):
        """Constructor.

        Parameters
        ----------
        species : list
            A list of strings, but not of ``Species``.
            The strings suggest serials of ``Species`` to be observed.

        """
        pass  # XXX: Only used for doc string

    def __cinit__(self, species):
        cdef vector[string] cpp_species
        for serial in species:
            cpp_species.push_back(tostring(serial))
        self.thisptr = new shared_ptr[Cpp_NumberObserver](
            new Cpp_NumberObserver(cpp_species))

    def __dealloc__(self):
        del self.thisptr

    def num_steps(self):
        """Return the number of steps."""
        return self.thisptr.get().num_steps()

    def next_time(self):
        """Return the next time for logging."""
        return self.thisptr.get().next_time()

    def data(self):
        """Return a list of the numbers of molecules you specified.

        Returns
        -------
        list:
            A list of lists of the number of molecules.
            The size of a return value is equal to ``num_steps``.
            Each element of a return value is a list consisting of
            time and the number of molecules specified at the construction.

        """
        cdef vector[vector[Real]] d = self.thisptr.get().data()
        retval = []
        cdef vector[vector[Real]].iterator it = d.begin()
        while it != d.end():
            retval.append(deref(it))
            inc(it)
        return retval

    def targets(self):
        """Return a list of ``Species``, which this ``Observer`` observes

        Returns
        -------
        list:
            A list of ``Species``. This is generated from arguments
            you gave at the construction.

        """
        cdef vector[Cpp_Species] species = self.thisptr.get().targets()

        retval = []
        cdef vector[Cpp_Species].iterator it = species.begin()
        while it != species.end():
            retval.append(
                 Species_from_Cpp_Species(
                     <Cpp_Species*>(address(deref(it)))))
            inc(it)
        return retval

    def save(self, filename):
        """Save data to an output with the given filename."""
        self.thisptr.get().save(tostring(filename))

    def as_base(self):
        """Clone self as a base class. This function is for developers."""
        retval = Observer()
        del retval.thisptr
        retval.thisptr = new shared_ptr[Cpp_Observer](
            <shared_ptr[Cpp_Observer]>deref(self.thisptr))
        return retval

    def reset(self):
        """Reset the internal state."""
        self.thisptr.get().reset()

cdef class TimingNumberObserver:
    """An ``Observer``class to log the number of molecules just at the time
    you assigned.

    TimingNumberObserver(t, species)

    """

    def __init__(self, vector[double] t, species):  #XXX: vector[Real]
        """Constructor.

        Parameters
        ----------
        t : list
            A list of times for logging. A time prior to the current
            time will be ignored.
        species : list
            A list of strings, but not of ``Species``.
            The strings suggest serials of ``Species`` to be observed.

        """
        pass  # XXX: Only used for doc string

    def __cinit__(self, vector[double] t, species):  #XXX: vector[Real]
        cdef vector[string] cpp_species
        for serial in species:
            cpp_species.push_back(tostring(serial))
        self.thisptr = new shared_ptr[Cpp_TimingNumberObserver](
            new Cpp_TimingNumberObserver(t, cpp_species))

    def __dealloc__(self):
        del self.thisptr

    def next_time(self):
        """Return the next time for logging."""
        return self.thisptr.get().next_time()

    def num_steps(self):
        """Return the number of steps."""
        return self.thisptr.get().num_steps()

    def data(self):
        """Return a list of the numbers of molecules you specified.

        Returns
        -------
        list:
            A list of lists of the number of molecules.
            The size of a return value is equal to ``num_steps``.
            Each element of a return value is a list consisting of
            time and the number of molecules specified at the construction.

        """
        cdef vector[vector[Real]] d = self.thisptr.get().data()
        retval = []
        cdef vector[vector[Real]].iterator it = d.begin()
        while it != d.end():
            retval.append(deref(it))
            inc(it)
        return retval

    def targets(self):
        """Return a list of ``Species``, which this ``Observer`` observes

        Returns
        -------
        list:
            A list of ``Species``. This is generated from arguments
            you gave at the construction.

        """
        cdef vector[Cpp_Species] species = self.thisptr.get().targets()

        retval = []
        cdef vector[Cpp_Species].iterator it = species.begin()
        while it != species.end():
            retval.append(
                 Species_from_Cpp_Species(
                     <Cpp_Species*>(address(deref(it)))))
            inc(it)
        return retval

    def save(self, filename):
        """Save data to an output with the given filename."""
        self.thisptr.get().save(tostring(filename))

    def as_base(self):
        """Clone self as a base class. This function is for developers."""
        retval = Observer()
        del retval.thisptr
        retval.thisptr = new shared_ptr[Cpp_Observer](
            <shared_ptr[Cpp_Observer]>deref(self.thisptr))
        return retval

    def reset(self):
        """Reset the internal state."""
        self.thisptr.get().reset()

cdef class FixedIntervalHDF5Observer:
    """An ``Observer`` class to log the state of ``World`` in HDF5 format
    with the fixed step interval.
    This ``Observer`` saves the ``World`` at the current time first, and
    then keeps saving every after the interval.

    FixedIntervalHDF5Observer(dt, filename)

    """

    def __init__(self, Real dt, filename):
        """Constructor.

        Parameters
        ----------
        dt : float
            A step interval for logging.
        filename : str
            A file name to be saved. Data are saved in HDF5 format.
            The extension name is recommended to be `.h5`.
            The file name can contain at most one formatting string like
            `%02d`, which will be replaced with the number of steps.
            When the file name contains no formmating string, data will
            be overwritten in a single file at every steps.

        """
        pass  # XXX: Only used for doc string

    def __cinit__(self, Real dt, filename):
        self.thisptr = new shared_ptr[Cpp_FixedIntervalHDF5Observer](
            new Cpp_FixedIntervalHDF5Observer(dt, tostring(filename)))

    def __dealloc__(self):
        del self.thisptr

    def next_time(self):
        """Return the next time for logging."""
        return self.thisptr.get().next_time()

    def num_steps(self):
        """Return the number of steps."""
        return self.thisptr.get().num_steps()

    def prefix(self):
        """Return a prefix of a file name given at the construction"""
        return self.thisptr.get().prefix().decode('UTF-8')

    def filename(self, idx=None):
        """Return a file name to be saved at the next time"""
        if idx is None:
            return self.thisptr.get().filename().decode('UTF-8')
        else:
            return self.thisptr.get().filename(<Integer>idx).decode('UTF-8')

    def as_base(self):
        """Clone self as a base class. This function is for developers."""
        retval = Observer()
        del retval.thisptr
        retval.thisptr = new shared_ptr[Cpp_Observer](
            <shared_ptr[Cpp_Observer]>deref(self.thisptr))
        return retval

    def reset(self):
        """Reset the internal state."""
        self.thisptr.get().reset()

cdef class FixedIntervalCSVObserver:
    """An ``Observer`` class to log the state of ``World`` in CSV format
    with the fixed step interval.
    This ``Observer`` saves the ``World`` at the current time first, and
    then keeps saving every after the interval.

    FixedIntervalCSVObserver(dt, filename, species=None)

    """

    def __init__(self, Real dt, filename, species=None):
        """Constructor.

        Parameters
        ----------
        dt : float
            A step interval for logging.
        filename : str
            A file name to be saved. Data are saved in CSV format.
            The extension name is recommended to be `.csv` or `.txt`.
            The file name can contain at most one formatting string like
            `%02d`, which will be replaced with the number of steps.
            When the file name contains no formmating string, data will
            be overwritten in a single file at every steps.
            The first line in a file represents labels for each row.
            Each column contains a position, a radius, and a serial id
            for the ``Species``.
        species : list
            A list of strings, but not of ``Species``.
            The strings suggest serials of ``Species`` to be observed.

        """
        pass  # XXX: Only used for doc string

    def __cinit__(self, Real dt, filename, species=None):
        cdef vector[string] cpp_species
        if species is None:
            self.thisptr = new shared_ptr[Cpp_FixedIntervalCSVObserver](
                new Cpp_FixedIntervalCSVObserver(dt, tostring(filename)))
        else:
            for serial in species:
                cpp_species.push_back(tostring(serial))
            self.thisptr = new shared_ptr[Cpp_FixedIntervalCSVObserver](
                new Cpp_FixedIntervalCSVObserver(
                    dt, tostring(filename), cpp_species))

    def __dealloc__(self):
        del self.thisptr

    def next_time(self):
        """Return the next time for logging."""
        return self.thisptr.get().next_time()

    def num_steps(self):
        """Return the number of steps."""
        return self.thisptr.get().num_steps()

    def log(self, w):
        """Force to log the given ``World`` to a file.

        Parameters
        ----------
        w : Space
            A ``Space`` (``World``) to be logged.

        Examples
        --------
        This is an easy way to save a ``World`` in CSV format without
        running a simulation.

        >>> w = spatiocyte.SpatiocyteWorld(Real3(1, 1, 1), 0.005)
        >>> w.bind_to(NetworkModel())
        >>> w.add_molecules(Species("A"), 3)
        >>> FixedIntervalCSVObserver(1, "test.csv").log(w)
        >>> print(open("test.csv").read())
        x,y,z,r,sid
        0.10614455552060439,0.66106605822212161,0.81500000000000006,0.0050000000000000001,0
        0.38375339303603129,0.37527767497325676,0.23999999999999999,0.0050000000000000001,0
        0.25311394008759508,0.05484827557301445,0.495,0.0050000000000000001,0
        """
        cdef Space space = w.as_base()
        # self.thisptr.get().log(space.thisptr.get())
        self.thisptr.get().log(deref(space.thisptr))

    def filename(self):
        """Return a file name to be saved at the next time"""
        return self.thisptr.get().filename().decode('UTF-8')

    def as_base(self):
        """Clone self as a base class. This function is for developers."""
        retval = Observer()
        del retval.thisptr
        retval.thisptr = new shared_ptr[Cpp_Observer](
            <shared_ptr[Cpp_Observer]>deref(self.thisptr))
        return retval

    def reset(self):
        """Reset the internal state."""
        self.thisptr.get().reset()

    def set_header(self, header):
        """Set the header. 'x,y,z,r,sid' as a default."""
        self.thisptr.get().set_header(tostring(header))

    def set_formatter(self, formatter):
        """Set the formatter.

        8 arguments are given: (1) the current time, (2-4) x, y, z,
        (5) radius, (6-7) lot and serial of a ParticleID, (8) species index.
        (1-5) are Real, and (6-8) Integer.
        Its default value is '%2%,%3%,%4%,%5%,%8%'.

        """
        self.thisptr.get().set_formatter(tostring(formatter))

cdef class CSVObserver:
    """An ``Observer`` class to log the state of ``World`` in CSV format.
    This ``Observer`` saves the ``World`` at the current time first, and
    then keeps saving every after steps.

    CSVObserver(filename, species=None)

    """

    def __init__(self, filename, species=None):
        """Constructor.

        Parameters
        ----------
        filename : str
            A file name to be saved. Data are saved in CSV format.
            The extension name is recommended to be `.csv` or `.txt`.
            The file name can contain at most one formatting string like
            `%02d`, which will be replaced with the number of steps.
            When the file name contains no formmating string, data will
            be overwritten in a single file at every steps.
            The first line in a file represents labels for each row.
            Each column contains a position, a radius, and a serial id
            for the ``Species``.
        species : list
            A list of strings, but not of ``Species``.
            The strings suggest serials of ``Species`` to be observed.

        """
        pass  # XXX: Only used for doc string

    def __cinit__(self, filename, species=None):
        cdef vector[string] cpp_species
        if species is None:
            self.thisptr = new shared_ptr[Cpp_CSVObserver](
                new Cpp_CSVObserver(tostring(filename)))
        else:
            for serial in species:
                cpp_species.push_back(tostring(serial))
            self.thisptr = new shared_ptr[Cpp_CSVObserver](
                new Cpp_CSVObserver(
                    tostring(filename), cpp_species))

    def __dealloc__(self):
        del self.thisptr

    def next_time(self):
        """Return the next time for logging."""
        return self.thisptr.get().next_time()

    def num_steps(self):
        """Return the number of steps."""
        return self.thisptr.get().num_steps()

    def log(self, w):
        """Force to log the given ``World`` to a file.

        Parameters
        ----------
        w : Space
            A ``Space`` (``World``) to be logged.

        Examples
        --------
        This is an easy way to save a ``World`` in CSV format without
        running a simulation.

        >>> w = spatiocyte.SpatiocyteWorld(Real3(1, 1, 1), 0.005)
        >>> w.bind_to(NetworkModel())
        >>> w.add_molecules(Species("A"), 3)
        >>> CSVObserver(1, "test.csv").log(w)
        >>> print(open("test.csv").read())
        x,y,z,r,sid
        0.10614455552060439,0.66106605822212161,0.81500000000000006,0.0050000000000000001,0
        0.38375339303603129,0.37527767497325676,0.23999999999999999,0.0050000000000000001,0
        0.25311394008759508,0.05484827557301445,0.495,0.0050000000000000001,0
        """
        cdef Space space = w.as_base()
        # self.thisptr.get().log(space.thisptr.get())
        self.thisptr.get().log(deref(space.thisptr))

    def filename(self):
        """Return a file name to be saved at the next time"""
        return self.thisptr.get().filename().decode('UTF-8')

    def as_base(self):
        """Clone self as a base class. This function is for developers."""
        retval = Observer()
        del retval.thisptr
        retval.thisptr = new shared_ptr[Cpp_Observer](
            <shared_ptr[Cpp_Observer]>deref(self.thisptr))
        return retval

    def reset(self):
        """Reset the internal state."""
        self.thisptr.get().reset()

    def set_header(self, header):
        """Set the header. 'x,y,z,r,sid' as a default."""
        self.thisptr.get().set_header(tostring(header))

    def set_formatter(self, formatter):
        """Set the formatter.

        8 arguments are given: (1) the current time, (2-4) x, y, z,
        (5) radius, (6-7) lot and serial of a ParticleID, (8) species index.
        (1-5) are Real, and (6-8) Integer.
        Its default value is '%2%,%3%,%4%,%5%,%8%'.

        """
        self.thisptr.get().set_formatter(tostring(formatter))

cdef class FixedIntervalTrajectoryObserver:
    """An ``Observer`` class to trace and log trajectories of diffusing
    particles in a ``World`` with the fixed step interval.
    This ``Observer`` logs at the current time first, and then keeps logging
    every after the interval.

    FixedIntervalTrajectoryObserver(dt, pids=None, resolve_boundary=None, subdt=None)

    """

    def __init__(self, Real dt, pids=None, resolve_boundary=None, subdt=None):
        """Constructor.

        Parameters
        ----------
        dt : float
            A step interval for logging.
        pids : list, optional
            A list of ``ParticleID``s.
        resolve_boundary : bool, optional
            If True, this ``Observer`` automatically resolves the effect
            of periodic boundary contidions by keeping shifts for each particles.
            Otherwise, this just logs positions within the size of ``World``
            with no care about boundary conditions.
        subdt : float, optional
            A step interval to check the periodic boundary.
            If None, use dt. This only has meaning when resolve_boundary is True.

        """
        pass  # XXX: Only used for doc string

    def __cinit__(self, Real dt, pids=None, resolve_boundary=None, subdt=None):
        cdef vector[Cpp_ParticleID] tmp

        if pids is None:
            self.thisptr = new shared_ptr[Cpp_FixedIntervalTrajectoryObserver](
                new Cpp_FixedIntervalTrajectoryObserver(
                    dt,
                    Cpp_FixedIntervalTrajectoryObserver.default_resolve_boundary() if resolve_boundary is None else <bool>resolve_boundary,
                    Cpp_FixedIntervalTrajectoryObserver.default_subdt() if subdt is None else <Real>subdt))
        else:
            for pid in pids:
                tmp.push_back(deref((<ParticleID>pid).thisptr))

            self.thisptr = new shared_ptr[Cpp_FixedIntervalTrajectoryObserver](
                new Cpp_FixedIntervalTrajectoryObserver(
                    dt, tmp,
                    Cpp_FixedIntervalTrajectoryObserver.default_resolve_boundary() if resolve_boundary is None else <bool>resolve_boundary,
                    Cpp_FixedIntervalTrajectoryObserver.default_subdt() if subdt is None else <Real>subdt))

    def __dealloc__(self):
        del self.thisptr

    def next_time(self):
        """Return the next time for logging."""
        return self.thisptr.get().next_time()

    def num_steps(self):
        """Return the number of steps."""
        return self.thisptr.get().num_steps()

    def num_tracers(self):
        """Return the number of tracer molecules."""
        return self.thisptr.get().num_tracers()

    def t(self):
        """Return time points at logging as a list."""
        return self.thisptr.get().t()

    def data(self):
        """Return a list of trajectories for each particles.

        Returns
        -------
        list:
            A list of lists of ``Real3``. An element of a return value
            is corresponding the trajectory of each particle. Thus, the size
            of a return value is the same with that of ``pids`` you gave
            at the construction.
            If a particle corresponding to the given ``ParticleID`` is missing,
            i.e. for a reaction, this ``Observer`` just skips to log the
            position. Therefore, lengths of the trajectories can be diverse.

        """
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
        """Clone self as a base class. This function is for developers."""
        retval = Observer()
        del retval.thisptr
        retval.thisptr = new shared_ptr[Cpp_Observer](
            <shared_ptr[Cpp_Observer]>deref(self.thisptr))
        return retval

    def reset(self):
        """Reset the internal state."""
        self.thisptr.get().reset()

cdef class TimingTrajectoryObserver:
    """An ``Observer`` class to trace and log trajectories of diffusing
    particles in a ``World`` at the given logging times.

    TimingTrajectoryObserver(t, pids, resolve_boundary=None, subdt=None)

    """

    def __init__(self, vector[double] t, pids=None, resolve_boundary=None, subdt=None):  # vector[Real]
        """Constructor.

        Parameters
        ----------
        t : list
            A list of the logging times.
            Times prior to the current time are ignored.
        pids : list, optional
            A list of ``ParticleID``s.
        resolve_boundary : bool, optional
            If True, this ``Observer`` automatically resolves the effect
            of periodic boundary contidions by keeping shifts for each particles.
            Otherwise, this just logs positions within the size of ``World``
            with no care about boundary conditions.
        subdt : float, optional
            A step interval to check the periodic boundary.
            If None, use dt. This only has meaning when resolve_boundary is True.

        """
        pass  # XXX: Only used for doc string

    def __cinit__(self, vector[double] t, pids=None, resolve_boundary=None, subdt=None):
        cdef vector[Cpp_ParticleID] tmp

        if pids is None:
            self.thisptr = new shared_ptr[Cpp_TimingTrajectoryObserver](
                new Cpp_TimingTrajectoryObserver(
                    t,
                    Cpp_TimingTrajectoryObserver.default_resolve_boundary() if resolve_boundary is None else <bool>resolve_boundary,
                    Cpp_TimingTrajectoryObserver.default_subdt() if subdt is None else <Real>subdt))
        else:
            for pid in pids:
                tmp.push_back(deref((<ParticleID>pid).thisptr))

            self.thisptr = new shared_ptr[Cpp_TimingTrajectoryObserver](
                new Cpp_TimingTrajectoryObserver(
                    t, tmp,
                    Cpp_TimingTrajectoryObserver.default_resolve_boundary() if resolve_boundary is None else <bool>resolve_boundary,
                    Cpp_TimingTrajectoryObserver.default_subdt() if subdt is None else <Real>subdt))

    def __dealloc__(self):
        del self.thisptr

    def next_time(self):
        """Return the next time for logging."""
        return self.thisptr.get().next_time()

    def num_steps(self):
        """Return the number of steps."""
        return self.thisptr.get().num_steps()

    def num_tracers(self):
        """Return the number of tracer molecules."""
        return self.thisptr.get().num_tracers()

    def t(self):
        """Return time points at logging as a list."""
        return self.thisptr.get().t()

    def data(self):
        """Return a list of trajectories for each particles.

        Returns
        -------
        list:
            A list of lists of ``Real3``. An element of a return value
            is corresponding the trajectory of each particle. Thus, the size
            of a return value is the same with that of ``pids`` you gave
            at the construction.
            If a particle corresponding to the given ``ParticleID`` is missing,
            i.e. for a reaction, this ``Observer`` just skips to log the
            position. Therefore, lengths of the trajectories can be diverse.

        """
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
        """Clone self as a base class. This function is for developers."""
        retval = Observer()
        del retval.thisptr
        retval.thisptr = new shared_ptr[Cpp_Observer](
            <shared_ptr[Cpp_Observer]>deref(self.thisptr))
        return retval

    def reset(self):
        """Reset the internal state."""
        self.thisptr.get().reset()

cdef class FixedIntervalTrackingObserver:
    """An ``Observer`` class to trace and log trajectories of diffusing
    particles in a ``World`` with the fixed step interval.
    This ``Observer`` logs at the current time first, and then keeps logging
    every after the interval.

    FixedIntervalTrackingObserver(dt, species, resolve_boundary=None, subdt=None, threshold=None)

    """

    def __init__(self, Real dt, species, resolve_boundary=None, subdt=None, threshold=None):
        """Constructor.

        Parameters
        ----------
        dt : float
            A step interval for logging.
        species : list
            A list of ``Species``.
        resolve_boundary : bool, optional
            If True, this ``Observer`` automatically resolves the effect
            of periodic boundary contidions by keeping shifts for each particles.
            Otherwise, this just logs positions within the size of ``World``
            with no care about boundary conditions.
        subdt : float, optional
            A step interval to check the periodic boundary.
            If None, use dt. This only has meaning when resolve_boundary is True.
        threshold : float, optional
            A maximum length to assume two particles are the same.
            If None, no threshold is used.

        """
        pass  # XXX: Only used for doc string

    def __init__(self, Real dt, species, resolve_boundary=None, subdt=None, threshold=None):
        cdef vector[Cpp_Species] tmp

        for sp in species:
            tmp.push_back(deref((<Species>sp).thisptr))

        self.thisptr = new shared_ptr[Cpp_FixedIntervalTrackingObserver](
            new Cpp_FixedIntervalTrackingObserver(
                dt, tmp,
                Cpp_FixedIntervalTrackingObserver.default_resolve_boundary() if resolve_boundary is None else <bool>resolve_boundary,
                Cpp_FixedIntervalTrackingObserver.default_subdt() if subdt is None else <Real>subdt,
                Cpp_FixedIntervalTrackingObserver.default_threshold() if threshold is None else <Real>threshold))

    def __dealloc__(self):
        del self.thisptr

    def next_time(self):
        """Return the next time for logging."""
        return self.thisptr.get().next_time()

    def num_steps(self):
        """Return the number of steps."""
        return self.thisptr.get().num_steps()

    def num_tracers(self):
        """Return the number of tracer molecules."""
        return self.thisptr.get().num_tracers()

    def t(self):
        """Return time points at logging as a list."""
        return self.thisptr.get().t()

    def data(self):
        """Return a list of trajectories for each particles.

        Returns
        -------
        list:
            A list of lists of ``Real3``. An element of a return value
            is corresponding the trajectory of each particle. Thus, the size
            of a return value is the same with that of ``pids`` you gave
            at the construction.
            If a particle corresponding to the given ``ParticleID`` is missing,
            i.e. for a reaction, this ``Observer`` just skips to log the
            position. Therefore, lengths of the trajectories can be diverse.

        """
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
        """Clone self as a base class. This function is for developers."""
        retval = Observer()
        del retval.thisptr
        retval.thisptr = new shared_ptr[Cpp_Observer](
            <shared_ptr[Cpp_Observer]>deref(self.thisptr))
        return retval

    def reset(self):
        """Reset the internal state."""
        self.thisptr.get().reset()

cdef class TimeoutObserver:
    """An ``Observer``class to stop simulation at the given calculation time.

    TimeoutObserver(interval)

    """

    def __init__(self, interval=None):
        """Constructor.

        Parameters
        ----------
        interval : float
            timeout in seconds.

        """
        pass  # XXX: Only used for doc string

    def __cinit__(self, interval=None):
        if interval is None:
            self.thisptr = new shared_ptr[Cpp_TimeoutObserver](
                new Cpp_TimeoutObserver())
        else:
            self.thisptr = new shared_ptr[Cpp_TimeoutObserver](
                new Cpp_TimeoutObserver(<Real>interval))

    def __dealloc__(self):
        del self.thisptr

    def interval(self):
        """Return the timeout in seconds."""
        return self.thisptr.get().interval()

    def duration(self):
        """Return the last time to be called."""
        return self.thisptr.get().duration()

    def accumulation(self):
        """Return the accumulation time."""
        return self.thisptr.get().accumulation()

    def as_base(self):
        """Clone self as a base class. This function is for developers."""
        retval = Observer()
        del retval.thisptr
        retval.thisptr = new shared_ptr[Cpp_Observer](
            <shared_ptr[Cpp_Observer]>deref(self.thisptr))
        return retval

    def reset(self):
        """Reset the internal state."""
        self.thisptr.get().reset()
