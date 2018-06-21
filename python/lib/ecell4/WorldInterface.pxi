cdef class WorldInterface:
    """An abstract base class of all worlds. This is for developers.

    WorldInterface()

    """

    def __init__(self):
        """Constructor"""
        pass

    def __cinit__(self):
        self.thisptr = new shared_ptr[Cpp_WorldInterface](
            <Cpp_WorldInterface*>(new Cpp_CompartmentSpaceVectorImpl(
                Cpp_Real3(1, 1, 1)))) #XXX: DUMMY

    def __dealloc__(self):
        del self.thisptr

    # Real t()
    def t(self):
        """Return the time of the world."""
        return self.thisptr.get().t()

    # void set_t(Real&)
    def set_t(self, Real t):
        """set_t(t)

        Set the value of the time of the world.

        Parameters
        ----------
        t : Real
            The time of the world

        """
        self.thisptr.get().set_t(t)

    # void save(string&)
    def save(self, filename):
        """save(filename)

        Save the world to a file.

        Parameters
        ----------
        filename : str
            a filename to save to

        """
        self.thisptr.get().save(tostring(filename))

    # void load(string&)
    def load(self, filename):
        """load(filename)

        Load the world from a file.

        Parameters
        ----------
        filename : str
            a filename to load from

        """
        self.thisptr.get().load(tostring(filename))

    # Real volume()
    def volume(self):
        """Return the volume of the world."""
        return self.thisptr.get().volume()

    # # Integer num_species()

    # # bool has_species(Cpp_Species&)
    # def has_species(self, Species sp):
    #     """has_species(sp) -> bool
    #
    #     Check if the given species is in the space or not.
    #
    #     Args:
    #         sp (Species): A species to be found.
    #
    #     Returns:
    #         bool: True if the species in the space.
    #
    #     """
    #     return self.thisptr.get().has_species(deref(sp.thisptr))

    # Integer num_molecules(Cpp_Species&)
    def num_molecules(self, Species sp):
        """num_molecules(sp) -> Integer

        Return the number of molecules.

        Parameters
        ----------
        sp : Species
            A species whose molecules you count

        Returns
        -------
        Integer:
            The number of molecules (of a given species)

        """
        return self.thisptr.get().num_molecules(deref(sp.thisptr))

    # Integer num_molecules_exact(Cpp_Species&)
    def num_molecules_exact(self, Species sp):
        """num_molecules_exact(sp) -> Integer

        Return the number of molecules of a given species.

        Parameters
        ----------
        sp : Species
            A species whose molecules you count

        Returns
        -------
        Integer:
            The number of molecules of a given species

        """
        return self.thisptr.get().num_molecules_exact(deref(sp.thisptr))

    # Real get_value(Cpp_Species &sp)
    def get_value(self, Species sp):
        """get_value(sp) -> Real

        Return the value (number) corresponding the given Species.

        Parameters
        ----------
        sp : Species
            a species whose value you require

        Returns
        -------
        Real:
            the value

        """
        return self.thisptr.get().get_value(deref(sp.thisptr))

    # Real get_value_exact(Cpp_Species &sp)
    def get_value_exact(self, Species sp):
        """get_value_exact(sp) -> Real

        Return the value (number) corresponding the given Species.

        Parameters
        ----------
        sp : Species
            a species whose value you require

        Returns
        -------
        Real:
            the value

        """
        return self.thisptr.get().get_value_exact(deref(sp.thisptr))

    # Cpp_Real3 edge_lengths()
    def edge_lengths(self):
        """edge_lengths() -> Real3

        Return the edge lengths of the world.

        """
        cdef Cpp_Real3 lengths = self.thisptr.get().edge_lengths()
        return Real3_from_Cpp_Real3(address(lengths))

    # Cpp_Real3 actual_lengths()
    def actual_lengths(self):
        """actual_lengths() -> Real3

        Return the actual edge lengths of the world.
        Same as ``edge_lengths``.
        """
        cdef Cpp_Real3 lengths = self.thisptr.get().actual_lengths()
        return Real3_from_Cpp_Real3(address(lengths))

    # Integer num_particles()
    # Integer num_particles(Cpp_Species&)
    def num_particles(self, Species sp = None):
        """num_particles(sp=None) -> Integer

        Return the number of particles.

        Parameters
        ----------
        sp : Species, optional
            The species of particles to count
            If no species is given, return the total number of particles.

        Returns
        -------
        Integer:
            The number of particles (of the given species)

        """
        if sp is None:
            return self.thisptr.get().num_particles()
        else:
            return self.thisptr.get().num_particles(deref(sp.thisptr))

    # Integer num_particles_exact(Cpp_Species&)
    def num_particles_exact(self, Species sp):
        """num_particles_exact(sp) -> Integer

        Return the number of particles of a given species.

        Parameters
        ----------
        sp : Species
            The species of particles to count

        Returns
        -------
        Integer:
            The number of particles of a given species

        """
        return self.thisptr.get().num_particles_exact(deref(sp.thisptr))

    # bool has_particle(Cpp_ParticleID&)
    def has_particle(self, ParticleID pid):
        """has_particle(pid) -> bool

        Check if a particle associated with a given particle id exists.

        Parameters
        ----------
        pid : ParticleID
            A particle id to check

        Returns
        -------
        bool:
            If a particle exists, return True. Otherwise return False

        """
        return self.thisptr.get().has_particle(deref(pid.thisptr))

    # pair[Cpp_ParticleID, Cpp_Particle] get_particle(Cpp_ParticleID&)
    def get_particle(self, ParticleID pid):
        """get_particle(pid) -> (ParticleID, Particle)

        Return the particle associated a given ParticleID.

        Parameters
        ----------
        pid : ParticleID
            An id of the particle you want

        Returns
        -------
        tuple:
            A pair of ParticleID and Particle

        """
        cdef pair[Cpp_ParticleID, Cpp_Particle] \
            pid_particle_pair = self.thisptr.get().get_particle(deref(pid.thisptr))
        return (ParticleID_from_Cpp_ParticleID(address(pid_particle_pair.first)),
                Particle_from_Cpp_Particle(address(pid_particle_pair.second)))

    # vector[pair[Cpp_ParticleID, Cpp_Particle]] list_particles()
    # vector[pair[Cpp_ParticleID, Cpp_Particle]] list_particles(Cpp_Species&)
    def list_particles(self, Species sp = None):
        """list_particles(sp) -> [(ParticleID, Particle)]

        Return the list of particles.

        Parameters
        ----------
        sp : Species, optional
            The species of particles to list up
            If no species is given, return the whole list of particles.

        Returns
        -------
        list:
            The list of particles (of the given species)

        """
        cdef vector[pair[Cpp_ParticleID, Cpp_Particle]] particles
        if sp is None:
            particles = self.thisptr.get().list_particles()
        else:
            particles = self.thisptr.get().list_particles(deref(sp.thisptr))

        retval = []
        cdef vector[pair[Cpp_ParticleID, Cpp_Particle]].iterator \
            it = particles.begin()
        while it != particles.end():
            retval.append(
                (ParticleID_from_Cpp_ParticleID(
                     <Cpp_ParticleID*>(address(deref(it).first))),
                 Particle_from_Cpp_Particle(
                     <Cpp_Particle*>(address(deref(it).second)))))
            inc(it)
        return retval

    # vector[pair[Cpp_ParticleID, Cpp_Particle]] list_particles_exact(Cpp_Species&)
    def list_particles_exact(self, Species sp):
        """list_particles_exact(sp) -> [(ParticleID, Particle)]

        Return the list of particles of a given species.

        Parameters
        ----------
        sp : Species
            The species of particles to list up

        Returns
        -------
        list:
            The list of particles of a given species

        """
        cdef vector[pair[Cpp_ParticleID, Cpp_Particle]] particles
        particles = self.thisptr.get().list_particles(deref(sp.thisptr))

        retval = []
        cdef vector[pair[Cpp_ParticleID, Cpp_Particle]].iterator \
            it = particles.begin()
        while it != particles.end():
            retval.append(
                (ParticleID_from_Cpp_ParticleID(
                     <Cpp_ParticleID*>(address(deref(it).first))),
                 Particle_from_Cpp_Particle(
                     <Cpp_Particle*>(address(deref(it).second)))))
            inc(it)
        return retval

cdef WorldInterface WorldInterface_from_Cpp_WorldInterface(shared_ptr[Cpp_WorldInterface] world):
    r = WorldInterface()
    r.thisptr.swap(world)
    return r
