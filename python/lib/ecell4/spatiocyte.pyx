import collections
from cython.operator cimport dereference as deref, preincrement as inc
from cython cimport address
from libcpp.string cimport string
from libcpp.vector cimport vector

from ecell4.types cimport *
from ecell4.shared_ptr cimport shared_ptr
from ecell4.optional cimport optional
from ecell4.core cimport *

from deprecated import deprecated

## ReactionInfo
cdef class ReactionInfo:
    """A class stores detailed information about a reaction in spatiocyte.

    ReactionInfo(t, reactants, products)

    """

    def __init__(self, Real t, reactants, products):
        """Constructor.

        Parameters
        ----------
        t : Real
            A time when a reaction occurred
        reactants : [ReactionInfoItem]
            A list of reactants.
        products : [ReactionInfoItem]
            A list of products.
        """
        pass  #XXX: only used for doc string


    def __cinit__(self, Real t, reactants, products):
        cdef vector[CppReactionInfoItem] cpp_reactants
        cdef vector[CppReactionInfoItem] cpp_products

        for item in reactants:
            cpp_reactants.push_back(deref((<ReactionInfoItem>item).thisptr))

        for item in products:
            cpp_products.push_back(deref((<ReactionInfoItem>item).thisptr))

        self.thisptr = new CppReactionInfo(t, cpp_reactants, cpp_products)

    def __dealloc__(self):
        del self.thisptr

    def t(self):
        """Return time when the reaction occurred."""
        return self.thisptr.t()

    def reactants(self):
        """Return a list of reactants

        Returns
        -------
        [ReactionInfoItem]:
            A list of information of the reactants
        """
        cdef vector[CppReactionInfoItem] particles = self.thisptr.reactants()

        retval = []
        cdef vector[CppReactionInfoItem].iterator itr = particles.begin()
        while itr != particles.end():
            retval.append(wrap_reaction_info_item(deref(itr)))
            inc(itr)

        return retval

    def products(self):
        """Return a list of products

        Returns
        -------
        [ReactionInfoItem]:
            A list of information of the products
        """
        cdef vector[CppReactionInfoItem] particles = self.thisptr.products()

        retval = []
        cdef vector[CppReactionInfoItem].iterator itr = particles.begin()
        while itr != particles.end():
            retval.append(wrap_reaction_info_item(deref(itr)))
            inc(itr)

        return retval

    def __reduce__(self):
        return (ReactionInfo, (self.t(), self.reactants(), self.products()))

cdef ReactionInfo ReactionInfo_from_Cpp_ReactionInfo(CppReactionInfo* ri):
    cdef CppReactionInfo *new_obj = new CppReactionInfo(<CppReactionInfo> deref(ri))
    r = ReactionInfo(0, [], [])
    del r.thisptr
    r.thisptr = new_obj
    return r

cdef ReactionInfoItem wrap_reaction_info_item(CppReactionInfoItem item):
    retval = ReactionInfoItem()
    del retval.thisptr
    retval.thisptr = new CppReactionInfoItem(item.pid, item.species, item.voxel)
    return retval

## Voxel
cdef class Voxel:

    def __cinit__(self):
        pass

    def __dealloc__(self):
        if self.thisptr:
            del self.thisptr

    def list_neighbors(self):
        ls = []
        for i in range(self.thisptr.num_neighbors()):
            ls.append(wrap_voxel(self.thisptr.get_neighbor(i)))
        return ls

    def position(self):
        cdef Cpp_Real3 position
        if self.thisptr:
            position = self.thisptr.position()
            return Real3_from_Cpp_Real3(address(position))
        return None

cdef Voxel wrap_voxel(CppVoxel voxel):
    cdef Voxel retval = Voxel()
    retval.thisptr = new CppVoxel(voxel)
    return retval

## SpatiocyteWorld
#  a python wrapper for Cpp_SpatiocyteWorld
cdef class SpatiocyteWorld:
    """A class containing the properties of the spatiocyte world.

    SpatiocyteWorld(edge_lengths=None, voxel_radius=None, GSLRandomNumberGenerator rng=None)

    """

    def __init__(self, edge_lengths = None, voxel_radius = None,
                 GSLRandomNumberGenerator rng = None):
        """Constructor.

        Parameters
        ----------
        edge_lengths : Real3, optional
            A size of the World.
        voxel_radius : Real, optional
            A radius of a voxel.
        rng : GSLRandomNumberGenerator, optional
            A random number generator.

        """
        pass

    def __cinit__(self, edge_lengths = None, voxel_radius = None,
                  GSLRandomNumberGenerator rng = None):
        cdef string filename

        if edge_lengths is None:
            self.thisptr = new shared_ptr[Cpp_SpatiocyteWorld](new Cpp_SpatiocyteWorld())
        elif voxel_radius is None:
            if isinstance(edge_lengths, Real3):
                self.thisptr = new shared_ptr[Cpp_SpatiocyteWorld](
                    new Cpp_SpatiocyteWorld(
                        deref((<Real3>edge_lengths).thisptr)))
            else:
                filename = tostring(edge_lengths)
                self.thisptr = new shared_ptr[Cpp_SpatiocyteWorld](
                    new Cpp_SpatiocyteWorld(filename))
        elif rng is None:
            self.thisptr = new shared_ptr[Cpp_SpatiocyteWorld](
                new Cpp_SpatiocyteWorld(
                    deref((<Real3>edge_lengths).thisptr), <Real>voxel_radius))
        else:
            self.thisptr = new shared_ptr[Cpp_SpatiocyteWorld](
                new Cpp_SpatiocyteWorld(
                    deref((<Real3>edge_lengths).thisptr), <Real>voxel_radius,
                    deref(rng.thisptr)))

    def __dealloc__(self):
        # XXX: Here, we release shared pointer,
        #      and if reference count to the SpatiocyteWorld object,
        #      it will be released automatically.
        del self.thisptr


    # WorldInterface

    def t(self):
        """Return the time of the world."""
        return self.thisptr.get().t()

    def set_t(self, Real t):
        """set_t(t)

        Set the value of the time of the world.

        Parameters
        ----------
        t : Real
            The time of the world

        """
        self.thisptr.get().set_t(t)

    def save(self, filename):
        """save(filename)

        Save the world to a file.

        Parameters
        ----------
        filename : str
            A filename to save to

        """
        self.thisptr.get().save(tostring(filename))

    def load(self, filename):
        """load(filename)

        Load the world from a file.

        Parameters
        ----------
        filename : str
            A filename to load from

        """
        self.thisptr.get().load(tostring(filename))

    def volume(self):
        """Return the volume of the world."""
        return self.thisptr.get().volume()

    def has_species(self, Species species):
        """has_species(species)

        Check if the world has a given species

        Parameters
        ----------
        species: Species
            A species to check existance in the world.

        Returns
        -------
        bool:
            if a species exists, this is true. Otherwise false
        """
        return self.thisptr.get().has_species(deref(species.thisptr))

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
        # if sp is None:
        #     return self.thisptr.get().num_molecules()
        # else:
        #     return self.thisptr.get().num_molecules(deref(sp.thisptr))
        return self.thisptr.get().num_molecules(deref(sp.thisptr))

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

    def edge_lengths(self):
        """edge_lengths() -> Real3

        Return the edge lengths of the world.

        """
        cdef Cpp_Real3 lengths = self.thisptr.get().edge_lengths()
        return Real3_from_Cpp_Real3(address(lengths))

    @deprecated(suggest="edge_lengths()")
    def actual_lengths(self):
        return self.edge_lengths()

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
            if a particle exists, this is true. Otherwise false

        """
        return self.thisptr.get().has_particle(deref(pid.thisptr))

    def get_particle(self, ParticleID pid):
        """get_particle(pid) -> (ParticleID, Particle)

        Return the particle associated a given ParticleID.

        Parameters
        ----------
        pid : ParticleID
            A id of the particle you want

        Returns
        -------
        tuple:
            A pair of ParticleID and Particle

        """
        cdef pair[Cpp_ParticleID, Cpp_Particle] \
            pid_particle_pair = self.thisptr.get().get_particle(deref(pid.thisptr))
        return (ParticleID_from_Cpp_ParticleID(address(pid_particle_pair.first)),
                Particle_from_Cpp_Particle(address(pid_particle_pair.second)))

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
        particles = self.thisptr.get().list_particles_exact(deref(sp.thisptr))

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


    # Particle Manipulation

    def new_particle(self, arg1, Real3 arg2=None):
        """new_particle(arg1, arg2=None) -> (ParticleID, Particle)

        Create a new particle.

        Parameters
        ----------
        arg1 : Particle
            A particle to be placed.

        or

        arg1 : Species
            A species of a particle
        arg2 : Real3
            A coordinate to place a particle

        Returns
        -------
        output: ParticleID or None
            A ParticleID of the new particle

        """
        cdef optional[Cpp_ParticleID] pid

        if arg2 is None:
            pid = self.thisptr.get().new_particle(deref((<Particle> arg1).thisptr))
        else:
            pid = self.thisptr.get().new_particle(deref((<Species> arg1).thisptr), deref(arg2.thisptr))
        if pid.is_initialized():
            return ParticleID_from_Cpp_ParticleID(address(pid.get()))

        return None

    def remove_particle(self, ParticleID pid):
        """remove_particle(pid)

        Remove the particle associated with a given ParticleID.

        Parameters
        ----------
        pid : ParticleID
            A id of particle to remove

        """
        self.thisptr.get().remove_particle(deref(pid.thisptr))

    def list_structure_particles(self):
        """list_strucutre_particles() -> [(ParticleID, Particle)]

        Return the list of structure particles

        Returns
        -------
        list:
            The list of particles constructing a structure
        """
        cdef vector[pair[Cpp_ParticleID, Cpp_Particle]] particles
        particles = self.thisptr.get().list_structure_particles()

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

    def list_non_structure_particles(self):
        """list_strucutre_particles() -> [(ParticleID, Particle)]

        Return the list of non-structure particles

        Returns
        -------
        list:
            The list of particles not constructing a structure
        """
        cdef vector[pair[Cpp_ParticleID, Cpp_Particle]] particles
        particles = self.thisptr.get().list_non_structure_particles()

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

    def update_particle(self, ParticleID pid, Particle p):
        """update_particle(pid, p)

        Update a particle.

        Parameters
        ----------
        pid : ParticleID
            A particle id of the particle to update
        p : Particle
            The information to update a particle

        Returns
        -------
        bool:
            True if a new particle was created.

        """
        return self.thisptr.get().update_particle(deref(pid.thisptr), deref(p.thisptr))


    # Molecule Manipulation

    def add_molecules(self, Species sp, Integer num, shape=None):
        """add_molecules(sp, num, shape=None)

        Add some molecules.

        Parameters
        ----------
        sp : Species
            A species of molecules to add
        num : Integer
            The number of molecules to add
        shape : Shape, optional
            A shape to add molecules on

        """
        if shape is None:
            self.thisptr.get().add_molecules(deref(sp.thisptr), num)
        else:
            self.thisptr.get().add_molecules(
                deref(sp.thisptr), num, deref((<Shape>(shape.as_base())).thisptr))

    def remove_molecules(self, Species sp, Integer num):
        """remove_molecules(sp, num)

        Remove the molecules.

        Parameters
        ----------
        sp : Species
            A species whose molecules to remove
        num : Integer
            A number of molecules to be removed

        """
        self.thisptr.get().remove_molecules(deref(sp.thisptr), num)


    # SpatiocyteWorld API

    @deprecated(suggest="volume()")
    def actual_volume(self):
        return self.volume()

    def voxel_volume(self):
        """Return the volume of a voxel."""
        return self.thisptr.get().voxel_volume()

    def get_volume(self, Species sp):
        """get_volume(sp) -> Real

        Return a volume of the given structure.

        Parameters
        ----------
        sp : Species
            A species for the target structure.

        Returns
        -------
        Real:
            A total volume of voxels belonging to the structure.

        """
        return self.thisptr.get().get_volume(deref(sp.thisptr))

    def get_voxel(self, ParticleID pid):
        """get_voxel(pid)

        Return a voxel occupied with the particle associated with the given ParticleID.

        Parameters
        ----------
        pid : ParticleID
            An id of the particle occupying the voxel

        Returns
        -------
        output: ParticleVoxel or None
        """
        cdef optional[Cpp_ParticleVoxel] voxel
        voxel = self.thisptr.get().find_voxel(deref(pid.thisptr))

        if voxel.is_initialized():
            return ParticleVoxel_from_Cpp_ParticleVoxel(address(voxel.get()))

        return None

    def get_voxel_at(self, Voxel voxel):
        """get_voxel_at(voxel)

        Return a voxel at the given coordinate.

        Parameters
        ----------
        voxel: Voxel
            A voxel coordinate

        Returns
        -------
        output: (ParticleID, Species)
        """

        pid_species_pair = self.thisptr.get().get_voxel_at(deref(voxel.thisptr))
        return (ParticleID_from_Cpp_ParticleID(address(pid_species_pair.first)),
                Species_from_Cpp_Species(address(pid_species_pair.second)))

    @deprecated(suggest="remove_particle(pid)")
    def remove_voxel(self, ParticleID pid):
        self.thisptr.get().remove_voxel(deref(pid.thisptr))

    def set_value(self, Species sp, Real value):
        """set_value(sp, value)

        Set the value of the given species.

        Parameters
        ----------
        sp : Species
            a species whose value you set
        value : Real
            a value set

        """
        self.thisptr.get().set_value(deref(sp.thisptr), value)

    def list_species(self):
        """list_species() -> [Species]

        Return the list of species.

        Returns
        -------
        list:
            The list of species
        """
        cdef vector[Cpp_Species] species = self.thisptr.get().list_species()
        retval = []
        cdef vector[Cpp_Species].iterator it = species.begin()
        while it != species.end():
            retval.append(
                 Species_from_Cpp_Species(
                     <Cpp_Species*>(address(deref(it)))))
            inc(it)
        return retval

    @deprecated
    def list_structure_species(self):
        cdef vector[Cpp_Species] species = self.thisptr.get().list_structure_species()
        retval = []
        cdef vector[Cpp_Species].iterator it = species.begin()
        while it != species.end():
            retval.append(
                 Species_from_Cpp_Species(
                     <Cpp_Species*>(address(deref(it)))))
            inc(it)
        return retval

    @deprecated
    def list_non_structure_species(self):
        cdef vector[Cpp_Species] species = self.thisptr.get().list_non_structure_species()
        retval = []
        cdef vector[Cpp_Species].iterator it = species.begin()
        while it != species.end():
            retval.append(
                 Species_from_Cpp_Species(
                     <Cpp_Species*>(address(deref(it)))))
            inc(it)
        return retval

    @deprecated(suggest="num_particles(sp)")
    def num_voxels(self, Species sp = None):
        if sp is None:
            return self.thisptr.get().num_voxels()
        else:
            return self.thisptr.get().num_voxels(deref(sp.thisptr))

    @deprecated(suggest="num_particles_exact(sp)")
    def num_voxels_exact(self, Species sp):
        return self.thisptr.get().num_voxels_exact(deref(sp.thisptr))

    @deprecated(suggest="voxel.list_neighbors()")
    def get_neighbor(self, Voxel voxel, nrand):
        return wrap_voxel(voxel.thisptr.get_neighbor(nrand))

    def new_voxel(self, Species sp, Voxel voxel):
        """new_voxel(sp, voxel)

        Create a particle.

        Parameters
        ----------
        sp : Species
            A Species to put
        voxel : Voxel
            A location to put on

        Returns
        -------
        output: ParticleID or None
        """

        cdef optional[Cpp_ParticleID] pid
        pid = self.thisptr.get().new_voxel(deref(sp.thisptr), deref(voxel.thisptr));

        if pid.is_initialized():
            return ParticleID_from_Cpp_ParticleID(address(pid.get()))

        return None

    def new_voxel_structure(self, Species species, Voxel voxel):
        """new_voxel_structure(species, voxel)

        Create a particle.

        Parameters
        ----------
        species : Species
            The Species of particles to create
        voxel : Voxel
            A voxel to place the structure

        Returns
        -------
        output: ParticleID or None

        """
        cdef optional[Cpp_ParticleID] pid
        pid = self.thisptr.get().new_voxel_structure(deref(species.thisptr), deref(voxel.thisptr))

        if pid.is_initialized():
            return ParticleID_from_Cpp_ParticleID(address(pid.get()))

        return None

    @deprecated
    def update_voxel(self, ParticleID pid, ParticleVoxel v):
        return self.thisptr.get().update_voxel(deref(pid.thisptr), deref(v.thisptr))

    @deprecated
    def list_voxels(self, Species sp = None):
        cdef vector[pair[Cpp_ParticleID, Cpp_ParticleVoxel]] voxels
        if sp is None:
            voxels = self.thisptr.get().list_voxels()
        else:
            voxels = self.thisptr.get().list_voxels(deref(sp.thisptr))

        retval = []
        cdef vector[pair[Cpp_ParticleID, Cpp_ParticleVoxel]].iterator \
            it = voxels.begin()
        while it != voxels.end():
            retval.append(
                (ParticleID_from_Cpp_ParticleID(
                     <Cpp_ParticleID*>(address(deref(it).first))),
                 ParticleVoxel_from_Cpp_ParticleVoxel(
                     <Cpp_ParticleVoxel*>(address(deref(it).second)))))
            inc(it)
        return retval

    @deprecated
    def list_voxels_exact(self, Species sp):
        cdef vector[pair[Cpp_ParticleID, Cpp_ParticleVoxel]] voxels
        voxels = self.thisptr.get().list_voxels_exact(deref(sp.thisptr))

        retval = []
        cdef vector[pair[Cpp_ParticleID, Cpp_ParticleVoxel]].iterator \
            it = voxels.begin()
        while it != voxels.end():
            retval.append(
                (ParticleID_from_Cpp_ParticleID(
                     <Cpp_ParticleID*>(address(deref(it).first))),
                 ParticleVoxel_from_Cpp_ParticleVoxel(
                     <Cpp_ParticleVoxel*>(address(deref(it).second)))))
            inc(it)
        return retval

    @deprecated(suggest="has_particle(pid)")
    def has_voxel(self, ParticleID pid):
        return self.thisptr.get().has_voxel(deref(pid.thisptr))

    def voxel_radius(self):
        """Return the voxel radius."""
        return self.thisptr.get().voxel_radius()

    def size(self):
        """Return the size of voxels."""
        return self.thisptr.get().size()

    def shape(self):
        """shape() -> Integer3

        Return the triplet of sizes of column, row and layer.

        """
        cdef Cpp_Integer3 sizes = self.thisptr.get().shape()
        return Integer3_from_Cpp_Integer3(address(sizes))


    def bind_to(self, m):
        """bind_to(m)

        Bind a model to the world

        Parameters
        ----------
        m : Model
            A model to bind

        """
        self.thisptr.get().bind_to(Cpp_Model_from_Model(m))

    @deprecated(suggest="voxel.position()")
    def coordinate2position(self, Voxel voxel):
        return voxel.position()

    @deprecated(suggest="get_voxel_near_by(pos)")
    def position2coordinate(self, Real3 pos):
        return self.get_voxel_near(pos)

    def get_voxel_near_by(self, Real3 pos):
        """get_voxel_near_by(pos) -> Voxel

        Get a voxel near by a given position.

        Parameters
        ----------
        pos : Real3
            A position

        Returns
        -------
        Voxel:
            A voxel

        """
        return wrap_voxel(self.thisptr.get().position2voxel(deref(pos.thisptr)))

    def add_structure(self, Species sp, shape):
        """add_structure(sp, shape)

        Add a structure.

        Parameters
        ----------
        sp : Species
            A species suggesting the shape.
        shape : Shape
            A shape of the structure.

        """
        return self.thisptr.get().add_structure(
            deref(sp.thisptr), deref((<Shape>(shape.as_base())).thisptr))

    def rng(self):
        """Return a random number generator object."""
        return GSLRandomNumberGenerator_from_Cpp_RandomNumberGenerator(
            self.thisptr.get().rng())

    @staticmethod
    def calculate_voxel_volume(voxel_radius):
        """Calculate a voxel volume from a voxel radius."""
        return Cpp_SpatiocyteWorld.calculate_voxel_volume(voxel_radius)

    @staticmethod
    def calculate_hcp_lengths(voxel_radius):
        """calculate_hcp_lengths(Real voxel_radius) -> Real3

        Calculate HCP lengths (HCP_L, HCP_X, HCP_Y) from a voxel radius.

        """
        cdef Cpp_Real3 lengths = Cpp_SpatiocyteWorld.calculate_hcp_lengths(voxel_radius)
        return Real3_from_Cpp_Real3(address(lengths))

    @staticmethod
    def calculate_shape(Real3 edge_lengths, voxel_radius):
        """calculate_shape(Real3 edge_lengths, Real voxel_radius) -> Integer3

        Calculate World shape.

        """
        cdef Cpp_Integer3 shape = Cpp_SpatiocyteWorld.calculate_shape(
            deref(edge_lengths.thisptr), voxel_radius)
        return Integer3_from_Cpp_Integer3(address(shape))

    @staticmethod
    def calculate_volume(Real3 edge_lengths, voxel_radius):
        """calculate_volume(Real3 edge_lengths, Real voxel_radius) -> Real

        Calculate World volume.

        """
        return Cpp_SpatiocyteWorld.calculate_volume(
            deref(edge_lengths.thisptr), voxel_radius)

    def as_base(self):
        """Return self as a base class. Only for developmental use."""
        retval = WorldInterface()
        del retval.thisptr
        retval.thisptr = new shared_ptr[Cpp_WorldInterface](
            <shared_ptr[Cpp_WorldInterface]>deref(self.thisptr))
        return retval

cdef SpatiocyteWorld SpatiocyteWorld_from_Cpp_SpatiocyteWorld(
    shared_ptr[Cpp_SpatiocyteWorld] w):
    r = SpatiocyteWorld(Real3(1, 1, 1))
    r.thisptr.swap(w)
    return r

def create_spatiocyte_world_cell_list_impl(
    edge_lengths, voxel_radius, matrix_sizes, rng):
    cdef shared_ptr[Cpp_SpatiocyteWorld]* w = new shared_ptr[Cpp_SpatiocyteWorld](
        create_spatiocyte_world_cell_list_impl_alias(
            deref((<Real3>edge_lengths).thisptr), <Real>voxel_radius,
            deref((<Integer3>matrix_sizes).thisptr),
            deref((<GSLRandomNumberGenerator>rng).thisptr)))
    return SpatiocyteWorld_from_Cpp_SpatiocyteWorld(deref(w))

def create_spatiocyte_world_vector_impl(edge_lengths, voxel_radius, rng):
    cdef shared_ptr[Cpp_SpatiocyteWorld]* w = new shared_ptr[Cpp_SpatiocyteWorld](
        create_spatiocyte_world_vector_impl_alias(
            deref((<Real3>edge_lengths).thisptr), <Real>voxel_radius,
            deref((<GSLRandomNumberGenerator>rng).thisptr)))
    return SpatiocyteWorld_from_Cpp_SpatiocyteWorld(deref(w))

def create_spatiocyte_world_square_offlattice_impl(edge_length, voxel_radius, rng):
    cdef shared_ptr[Cpp_SpatiocyteWorld]* world = new shared_ptr[Cpp_SpatiocyteWorld](
        allocate_spatiocyte_world_square_offlattice_impl(
            <Real>edge_length, <Real>voxel_radius,
            deref((<GSLRandomNumberGenerator>rng).thisptr)))
    return SpatiocyteWorld_from_Cpp_SpatiocyteWorld(deref(world))

## SpatiocyteSimulator
#  a python wrapper for Cpp_SpatiocyteSimulator
cdef class SpatiocyteSimulator:
    """ A class running the simulation with the spatiocyte algorithm.

    SpatiocyteSimulator(m, w)

    """

    def __init__(self, SpatiocyteWorld w, m=None):
        """SpatiocyteSimulator(m, w)
        SpatiocyteSimulator(w)

        Constructor.

        Parameters
        ----------
        w : SpatiocyteWorld
            A world
        m : Model, optional
            A model

        """
        pass

    def __cinit__(self, SpatiocyteWorld w, m=None):
        if m is None:
            self.thisptr = new Cpp_SpatiocyteSimulator(deref(w.thisptr))
        else:
            self.thisptr = new Cpp_SpatiocyteSimulator(
                 deref(w.thisptr), Cpp_Model_from_Model(m))

    def __dealloc__(self):
        del self.thisptr

    def num_steps(self):
        """Return the number of steps."""
        return self.thisptr.num_steps()

    def step(self, upto = None):
        """step(upto=None) -> bool

        Step the simulation.

        Parameters
        ----------
        upto : Real, optional
            The time which to step the simulation up to

        Returns
        -------
        bool:
            True if the simulation did not reach the given time.
            When upto is not given, nothing will be returned.

        """
        if upto is None:
            self.thisptr.step()
        else:
            return self.thisptr.step(upto)

    def t(self):
        """Return the time."""
        return self.thisptr.t()

    def set_t(self, Real t_new):
        """set_t(t)

        Set the current time.

        Parameters
        ----------
        t : Real
            A current time.

        """
        self.thisptr.set_t(t_new)

    def dt(self):
        """Return the step interval."""
        return self.thisptr.dt()

    def next_time(self):
        """Return the scheduled time for the next step."""
        return self.thisptr.next_time()

    def set_dt(self, Real dt):
        """set_dt(dt)

        Set a step interval.

        Parameters
        ----------
        dt : Real
            A step interval

        """
        self.thisptr.set_dt(dt)

    def initialize(self):
        """Initialize the simulator."""
        self.thisptr.initialize()

    def check_reaction(self):
        """Return if any reaction occurred at the last step, or not."""
        return self.thisptr.check_reaction()

    def last_reactions(self):
        """last_reactions() -> [(ReactionRule, ReactionInfo)]

        Return reactions occuring at the last step.

        Returns
        -------
        list:
            The list of reaction rules and infos.

        """
        cdef vector[pair[Cpp_ReactionRule, CppReactionInfo]] reactions = self.thisptr.last_reactions()
        cdef vector[pair[Cpp_ReactionRule, CppReactionInfo]].iterator it = reactions.begin()
        retval = []
        while it != reactions.end():
            retval.append((
                ReactionRule_from_Cpp_ReactionRule(
                    <Cpp_ReactionRule*>(address(deref(it).first))),
                ReactionInfo_from_Cpp_ReactionInfo(
                    <CppReactionInfo*>(address(deref(it).second)))))
            inc(it)
        return retval

    def model(self):
        """Return the model bound."""
        return Model_from_Cpp_Model(self.thisptr.model())

    def world(self):
        """Return the world bound."""
        return SpatiocyteWorld_from_Cpp_SpatiocyteWorld(self.thisptr.world())

    def run(self, Real duration, observers=None, is_dirty=None):
        """run(duration, observers, is_dirty)

        Run the simulation.

        Parameters
        ----------
        duration : Real
            A duration for running a simulation.
            A simulation is expected to be stopped at t() + duration.
        observers : list of Obeservers, optional
            observers
        is_dirty : bool, default True
            If True, call initialize before running.

        """
        cdef vector[shared_ptr[Cpp_Observer]] tmp

        if is_dirty is None:
            if observers is None:
                self.thisptr.run(duration)
            elif isinstance(observers, collections.Iterable):
                for obs in observers:
                    tmp.push_back(deref((<Observer>(obs.as_base())).thisptr))
                self.thisptr.run(duration, tmp)
            else:
                self.thisptr.run(duration,
                    deref((<Observer>(observers.as_base())).thisptr))
        else:
            if observers is None:
                self.thisptr.run(duration, is_dirty)
            elif isinstance(observers, collections.Iterable):
                for obs in observers:
                    tmp.push_back(deref((<Observer>(obs.as_base())).thisptr))
                self.thisptr.run(duration, tmp, is_dirty)
            else:
                self.thisptr.run(duration,
                    deref((<Observer>(observers.as_base())).thisptr), is_dirty)

cdef SpatiocyteSimulator SpatiocyteSimulator_from_Cpp_SpatiocyteSimulator(Cpp_SpatiocyteSimulator* s):
    r = SpatiocyteSimulator(
        Model_from_Cpp_Model(s.model()), SpatiocyteWorld_from_Cpp_SpatiocyteWorld(s.world()))
    del r.thisptr
    r.thisptr = s
    return r

## SpatiocyteFactory
#  a python wrapper for Cpp_SpatiocyteFactory
cdef class SpatiocyteFactory:
    """ A factory class creating a SpatiocyteWorld instance and a SpatiocyteSimulator instance.

    SpatiocyteFactory(Real voxel_radius)

    """

    def __init__(self, voxel_radius=None):
        """SpatiocyteFactory(Real voxel_radius=None)

        Constructor.

        Parameters
        ----------
        voxel_radius : Real, optional
            A radius of a voxel.

        """
        pass

    def __cinit__(self, voxel_radius=None):
        self.thisptr = new Cpp_SpatiocyteFactory(
            Cpp_SpatiocyteFactory.default_voxel_radius() if voxel_radius is None else <Real>voxel_radius)

    def __dealloc__(self):
        del self.thisptr

    def rng(self, GSLRandomNumberGenerator rng):
        """rng(GSLRandomNumberGenerator) -> SpatiocyteFactory

        Set a random number generator, and return self.

        """
        cdef Cpp_SpatiocyteFactory *ptr = self.thisptr.rng_ptr(deref(rng.thisptr))
        assert ptr == self.thisptr
        return self

    def world(self, arg1=None):
        """world(arg1=None) -> SpatiocyteWorld

        Return a SpatiocyteWorld instance.

        Parameters
        ----------
        arg1 : Real3
            The lengths of edges of a SpatiocyteWorld created

        or

        arg1 : str
            The path of a HDF5 file for SpatiocyteWorld

        Returns
        -------
        SpatiocyteWorld:
            The created world

        """
        if arg1 is None:
            return SpatiocyteWorld_from_Cpp_SpatiocyteWorld(
                shared_ptr[Cpp_SpatiocyteWorld](self.thisptr.world()))
        elif isinstance(arg1, Real3):
            return SpatiocyteWorld_from_Cpp_SpatiocyteWorld(
                shared_ptr[Cpp_SpatiocyteWorld](
                    self.thisptr.world(deref((<Real3>arg1).thisptr))))
        elif isinstance(arg1, str):
            return SpatiocyteWorld_from_Cpp_SpatiocyteWorld(
                shared_ptr[Cpp_SpatiocyteWorld](self.thisptr.world(<string>(arg1))))
        else:
            return SpatiocyteWorld_from_Cpp_SpatiocyteWorld(
                shared_ptr[Cpp_SpatiocyteWorld](self.thisptr.world(
                    Cpp_Model_from_Model(arg1))))

    def simulator(self, SpatiocyteWorld arg1, arg2=None):
        """simulator(arg1, arg2) -> SpatiocyteSimulator

        Return a SpatiocyteSimulator instance.

        Parameters
        ----------
        arg1 : SpatiocyteWorld
            A world
        arg2 : Model, optional
            A simulation model

        Returns
        -------
        SpatiocyteSimulator:
            The created simulator

        """
        if arg2 is None:
            return SpatiocyteSimulator_from_Cpp_SpatiocyteSimulator(
                self.thisptr.simulator(deref(arg1.thisptr)))
        else:
            return SpatiocyteSimulator_from_Cpp_SpatiocyteSimulator(
                self.thisptr.simulator(
                    deref(arg1.thisptr), Cpp_Model_from_Model(arg2)))

    def create_world(self, arg1=None):
        """create_world(arg1=None) -> SpatiocyteWorld

        Return a SpatiocyteWorld instance.

        Parameters
        ----------
        arg1 : Real3
            The lengths of edges of a SpatiocyteWorld created

        or

        arg1 : str
            The path of a HDF5 file for SpatiocyteWorld

        Returns
        -------
        SpatiocyteWorld:
            The created world

        """
        import warnings; warnings.warn("Function 'create_world()' has moved to 'world()'", DeprecationWarning)
        return self.world(arg1)

    def create_simulator(self, SpatiocyteWorld arg1, arg2=None):
        """create_simulator(arg1, arg2) -> SpatiocyteSimulator

        Return a SpatiocyteSimulator instance.

        Parameters
        ----------
        arg1 : SpatiocyteWorld
            A world
        arg2 : Model, optional
            A simulation model

        Returns
        -------
        SpatiocyteSimulator:
            The created simulator

        """
        import warnings; warnings.warn("Function 'create_simulator()' has moved to 'simulator()'", DeprecationWarning)
        return self.simulator(arg1, arg2)

Factory = SpatiocyteFactory  # This is an alias
