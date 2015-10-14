import collections
from cython.operator cimport dereference as deref, preincrement as inc
from cython cimport address
from libcpp.string cimport string
from libcpp.vector cimport vector

from ecell4.types cimport *
from ecell4.shared_ptr cimport shared_ptr
from ecell4.core cimport *

## ReactionInfo
cdef class ReactionInfo:
    """A class stores detailed information about a reaction in lattice.

    ReactionInfo(t, reactants, products)

    """

    def __init__(self, Real t, reactants, products):
        """Constructor.

        Args:
          t (Real): A time when a reaction occurred
          reactants (list): A list of reactants.
            Reactants are given as a pair of ``ParticleID`` and ``Voxel``.
          products (list): A list of products.
            Products are given as a pair of ``ParticleID`` and ``Voxel``.

        """
        pass  #XXX: only used for doc string


    def __cinit__(self, Real t, reactants, products):
        cdef vector[pair[Cpp_ParticleID, Cpp_Voxel]] reactants_
        cdef vector[pair[Cpp_ParticleID, Cpp_Voxel]] products_

        for pid, p in reactants:
            reactants_.push_back(
                pair[Cpp_ParticleID, Cpp_Voxel](
                    deref((<ParticleID>pid).thisptr), deref((<Voxel>p).thisptr)))
        for pid, p in products:
            products_.push_back(
                pair[Cpp_ParticleID, Cpp_Voxel](
                    deref((<ParticleID>pid).thisptr), deref((<Voxel>p).thisptr)))

        self.thisptr = new Cpp_ReactionInfo(t, reactants_, products_)

    def __dealloc__(self):
        del self.thisptr

    def t(self):
        """Return a time when a reaction occurred."""
        return self.thisptr.t()

    def reactants(self):
        """Return a list of reactants

        Returns:
            list: A list of pairs of ``ParticleID`` and ``Voxel``.

        """
        cdef vector[pair[Cpp_ParticleID, Cpp_Voxel]] particles
        particles = self.thisptr.reactants()

        retval = []
        cdef vector[pair[Cpp_ParticleID, Cpp_Voxel]].iterator \
            it = particles.begin()
        while it != particles.end():
            retval.append(
                (ParticleID_from_Cpp_ParticleID(
                     <Cpp_ParticleID*>(address(deref(it).first))),
                 Voxel_from_Cpp_Voxel(
                     <Cpp_Voxel*>(address(deref(it).second)))))
            inc(it)
        return retval

    def products(self):
        """Return a list of products

        Returns:
            list: A list of pairs of ``ParticleID`` and ``Voxel``.

        """
        cdef vector[pair[Cpp_ParticleID, Cpp_Voxel]] particles
        particles = self.thisptr.products()

        retval = []
        cdef vector[pair[Cpp_ParticleID, Cpp_Voxel]].iterator \
            it = particles.begin()
        while it != particles.end():
            retval.append(
                (ParticleID_from_Cpp_ParticleID(
                     <Cpp_ParticleID*>(address(deref(it).first))),
                 Voxel_from_Cpp_Voxel(
                     <Cpp_Voxel*>(address(deref(it).second)))))
            inc(it)
        return retval

cdef ReactionInfo ReactionInfo_from_Cpp_ReactionInfo(Cpp_ReactionInfo* ri):
    cdef Cpp_ReactionInfo *new_obj = new Cpp_ReactionInfo(<Cpp_ReactionInfo> deref(ri))
    r = ReactionInfo(0, [], [])
    del r.thisptr
    r.thisptr = new_obj
    return r

## LatticeWorld
#  a python wrapper for Cpp_LatticeWorld
cdef class LatticeWorld:
    """A class containing the properties of the lattice world.

    LatticeWorld(edge_lengths=None, voxel_radius=None, GSLRandomNumberGenerator rng=None)

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
            self.thisptr = new shared_ptr[Cpp_LatticeWorld](new Cpp_LatticeWorld())
        elif voxel_radius is None:
            if isinstance(edge_lengths, Real3):
                self.thisptr = new shared_ptr[Cpp_LatticeWorld](
                    new Cpp_LatticeWorld(
                        deref((<Real3>edge_lengths).thisptr)))
            else:
                filename = tostring(edge_lengths)
                self.thisptr = new shared_ptr[Cpp_LatticeWorld](
                    new Cpp_LatticeWorld(filename))
        elif rng is None:
            self.thisptr = new shared_ptr[Cpp_LatticeWorld](
                new Cpp_LatticeWorld(
                    deref((<Real3>edge_lengths).thisptr), <Real>voxel_radius))
        else:
            self.thisptr = new shared_ptr[Cpp_LatticeWorld](
                new Cpp_LatticeWorld(
                    deref((<Real3>edge_lengths).thisptr), <Real>voxel_radius,
                    deref(rng.thisptr)))

    def __dealloc__(self):
        # XXX: Here, we release shared pointer,
        #      and if reference count to the LatticeWorld object,
        #      it will be released automatically.
        del self.thisptr

    def set_t(self, Real t):
        """set_t(t)

        Set the value of the time of the world.

        Parameters
        ----------
        t : Real
            The time of the world

        """
        self.thisptr.get().set_t(t)

    def t(self):
        """Return the time of the world."""
        return self.thisptr.get().t()

    def volume(self):
        """Return the volume of the world."""
        return self.thisptr.get().volume()

    def voxel_volume(self):
        """Return the volume of a voxel."""
        return self.thisptr.get().voxel_volume()

    def get_volume(self):
        """Return the actual volume of the world."""
        return self.thisptr.get().get_volume()

    def actual_lengths(self):
        """Return the actual edge lengths of the world.

        Returns
        -------
        Real3:
            The actual edge lengths of the world

        """
        cdef Cpp_Real3 lengths = self.thisptr.get().actual_lengths()
        return Real3_from_Cpp_Real3(address(lengths))

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
        tuple:
            A pair of ParticleID and Particle of a new particle

        """
        cdef pair[pair[Cpp_ParticleID, Cpp_Particle], bool] retval

        if arg2 is None:
            retval = self.thisptr.get().new_particle(deref((<Particle> arg1).thisptr))
        else:
            retval = self.thisptr.get().new_particle(deref((<Species> arg1).thisptr), deref(arg2.thisptr))
        return ((ParticleID_from_Cpp_ParticleID(address(retval.first.first)), Particle_from_Cpp_Particle(address(retval.first.second))), retval.second)

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

    def get_voxel(self, arg):
        """get_voxel(arg) -> (ParticleID, Voxle)

        Return the voxel having a particle associated with a given ParticleID
        or coordinate.

        Parameters
        ----------
        arg : ParticleID or Integer
            An id or coordiante of the particle in the voxel you want

        Returns
        -------
        tuple:
            A pair of ParticleID and Voxel

        """
        cdef pair[Cpp_ParticleID, Cpp_Voxel] pid_voxel_pair
        if isinstance(arg, ParticleID):
            pid_voxel_pair = self.thisptr.get().get_voxel(deref((<ParticleID>arg).thisptr))
        else:
            pid_voxel_pair = self.thisptr.get().get_voxel(<Integer>arg)
        return (ParticleID_from_Cpp_ParticleID(address(pid_voxel_pair.first)),
                Voxel_from_Cpp_Voxel(address(pid_voxel_pair.second)))

    # def on_structure(self, Voxel v):
    #     """Check if the given voxel would be on the proper structure at the coordinate
    #     Args:
    #         v (Voxel): a voxel scheduled to be placed
    #     Returns:
    #         bool: if it is on the proper structure, or not
    #     """
    #     return self.thisptr.get().on_structure(deref((<Voxel>v).thisptr))

    def on_structure(self, Species sp, Integer coord):
        """on_structure(sp, coord) -> bool

        Check if the given species would be on the proper structure at the coordinate.

        Parameters
        ----------
        sp : Species
            A species scheduled to be placed
        coord : Integer
            A coordinate to be occupied

        Returns
        -------
        bool:
            if it is on the proper structure, or not

        """
        return self.thisptr.get().on_structure(deref(sp.thisptr), coord)

    def remove_particle(self, ParticleID pid):
        """remove_particle(pid)

        Remove the particle associated with a given ParticleID.

        Parameters
        ----------
        pid : ParticleID
            A id of particle to remove

        """
        self.thisptr.get().remove_particle(deref(pid.thisptr))

    def remove_voxel(self, ParticleID pid):
        """remove_voxel(pid)

        Remove the particle associated with a given ParticleID.

        Parameters
        ----------
        pid : ParticleID
            A id of particle to remove

        """
        self.thisptr.get().remove_voxel(deref(pid.thisptr))

    def edge_lengths(self):
        """edge_lengths() -> Real3

        Return the edge lengths of the world.

        """
        cdef Cpp_Real3 lengths = self.thisptr.get().edge_lengths()
        return Real3_from_Cpp_Real3(address(lengths))

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

    def num_voxels(self, Species sp = None):
        """num_voxels(sp=None) -> Integer

        Return the number of voxels.

        Parameters
        ----------
        sp : Species, optional
            The species of particles to count

        Returns
        -------
        Integer:
            The number of voxels (of the given species)

        """
        if sp is None:
            return self.thisptr.get().num_voxels()
        else:
            return self.thisptr.get().num_voxels(deref(sp.thisptr))

    def num_voxels_exact(self, Species sp):
        """num_voxels_exact(sp) -> Integer

        Return the number of voxels of a given species.

        Parameters
        ----------
        sp : Species
            The species of particles to count

        Returns
        -------
        Integer:
            The number of voxels of a given species

        """
        return self.thisptr.get().num_voxels_exact(deref(sp.thisptr))

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

    def get_neighbor(self, coord, nrand):
        """get_neighbor(coord, nrand) -> Integer

        Return the neighbor coordinate of a given coordinate.

        Parameters
        ----------
        coord : Integer
            A coordinate of a voxel
        nrand : Integer
            A key in the range from 0 to 11 to assign a neighbor voxel

        Returns
        -------
        Integer:
            The coordinate of the neighbor voxel

        """
        return self.thisptr.get().get_neighbor(coord, nrand)

    def get_neighbor_private(self, coord, nrand):
        """get_neighbor_private(coord, nrand) -> Integer

        Return the neighbor coordinate of a given coordinate in private.

        Parameters
        ----------
        coord : Integer
            A private coordinate of a voxel
        nrand : Integer
            A key in the range from 0 to 11 to assign a neighbor voxel

        Returns
        -------
        Integer:
            The private coordinate of the neighbor voxel

        """
        return self.thisptr.get().get_neighbor_private(coord, nrand)

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

    def new_voxel(self, arg1, arg2=None):
        """new_voxel(arg1, arg2) -> (ParticleID, Voxel)

        Create a particle.

        Parameters
        ----------
        arg1 : Voxel
            The information to create

        or

        arg1 : Species
            The Species of particles to create
        arg2 : Integer
            The number of particles(voxels)

        Returns
        -------
        tuple:
            A pair of ParticleID and Voxel

        """
        cdef pair[pair[Cpp_ParticleID, Cpp_Voxel], bool] retval

        if arg2 is None:
            retval = self.thisptr.get().new_voxel(deref((<Voxel> arg1).thisptr))
        else:
            retval = self.thisptr.get().new_voxel(deref((<Species> arg1).thisptr), <Integer> arg2)
        return ((ParticleID_from_Cpp_ParticleID(address(retval.first.first)), Voxel_from_Cpp_Voxel(address(retval.first.second))), retval.second)

    def update_voxel(self, ParticleID pid, Voxel v):
        """update_voxel(pid, v) -> bool

        Update a particle.

        Parameters
        ----------
        pid : ParticleID
            A particle id of the particle to update
        v : Voxel
            The information to update

        Returns
        -------
        bool:
            whether to succeed to update the particle

        """
        return self.thisptr.get().update_voxel(deref(pid.thisptr), deref(v.thisptr))

    def list_voxels(self, Species sp = None):
        """list_voxels(sp=None) -> [ParitcleID, Voxel]

        Returns the list of voxels.

        Parameters
        ----------
        sp : Species, optional
            A species of particles to list up.
            If no species is given, return a list of all voxels.

        Returns
        -------
        list:
            The list of the pair of ParticleID and Voxel

        """
        cdef vector[pair[Cpp_ParticleID, Cpp_Voxel]] voxels
        if sp is None:
            voxels = self.thisptr.get().list_voxels()
        else:
            voxels = self.thisptr.get().list_voxels(deref(sp.thisptr))

        retval = []
        cdef vector[pair[Cpp_ParticleID, Cpp_Voxel]].iterator \
            it = voxels.begin()
        while it != voxels.end():
            retval.append(
                (ParticleID_from_Cpp_ParticleID(
                     <Cpp_ParticleID*>(address(deref(it).first))),
                 Voxel_from_Cpp_Voxel(
                     <Cpp_Voxel*>(address(deref(it).second)))))
            inc(it)
        return retval

    def list_voxels_exact(self, Species sp):
        """list_voxels_exact(sp) -> [ParitcleID, Voxel]

        Returns the list of voxels.

        Parameters
        ----------
        sp : Species, optional
            A species of particles to list up.
            If no species is given, return a list of all voxels.

        Returns
        -------
        list:
            The list of the pair of ParticleID and Voxel

        """
        cdef vector[pair[Cpp_ParticleID, Cpp_Voxel]] voxels
        voxels = self.thisptr.get().list_voxels_exact(deref(sp.thisptr))

        retval = []
        cdef vector[pair[Cpp_ParticleID, Cpp_Voxel]].iterator \
            it = voxels.begin()
        while it != voxels.end():
            retval.append(
                (ParticleID_from_Cpp_ParticleID(
                     <Cpp_ParticleID*>(address(deref(it).first))),
                 Voxel_from_Cpp_Voxel(
                     <Cpp_Voxel*>(address(deref(it).second)))))
            inc(it)
        return retval

    def has_voxel(self, ParticleID pid):
        """has_voxel(pid) -> bool

        Check if a particle exists.

        Parameters
        ----------
        pid : ParticleID
            A particle id of the particle to check

        Returns
        -------
        bool:
            whether a particle associated with a given particle id exists

        """
        return self.thisptr.get().has_voxel(deref(pid.thisptr))

    def voxel_radius(self):
        """Return the voxel radius."""
        return self.thisptr.get().voxel_radius()

    def col_size(self):
        """Return the size of the column of the world."""
        return self.thisptr.get().col_size()

    def row_size(self):
        """Return the size of row of the world."""
        return self.thisptr.get().row_size()

    def layer_size(self):
        """Return the size of layer of the world."""
        return self.thisptr.get().layer_size()

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
        self.thisptr.get().bind_to(deref(Cpp_Model_from_Model(m)))

    def private2position(self, Integer coord):
        """private2position(coord) -> Real3

        Transform a private coordinate to a position.

        """
        cdef Cpp_Real3 pos = self.thisptr.get().private2position(coord)
        return Real3_from_Cpp_Real3(address(pos))

    def coordinate2position(self, Integer coord):
        """coordinate2position(coord) -> Real3

        Transform a coordinate to a position.

        """
        cdef Cpp_Real3 pos = self.thisptr.get().coordinate2position(coord)
        return Real3_from_Cpp_Real3(address(pos))

    def position2coordinate(self, Real3 pos):
        """position2coordinate(pos) -> Integer

        Transform a position to a coordinate.

        Parameters
        ----------
        pos : Real3
            A position

        Returns
        -------
        Integer:
            A coordinate

        """
        return self.thisptr.get().position2coordinate(
            deref(pos.thisptr))

    def private2coord(self, Integer coord):
        """Transform a private coordinate to a coordinate."""
        return self.thisptr.get().private2coord(coord)

    def coord2private(self, Integer coord):
        """Transform a coordinate to a private coordinate."""
        return self.thisptr.get().coord2private(coord)

    def global2coord(self, Integer3 g):
        """global2coord(g) -> Integer

        Transform a global coordinate to a coordinate.

        Parameters
        ----------
        g : Integer3
            A global coordinate

        Returns
        -------
        Integer:
            A coordinate

        """
        return self.thisptr.get().global2coord(deref(g.thisptr))

    def coord2global(self, Integer coord):
        """coord2global(coord) -> Integer3

        Transform a coordinate to a global coordinate.

        """
        cdef Cpp_Integer3 g = self.thisptr.get().coord2global(coord)
        return Integer3_from_Cpp_Integer3(address(g))

    def global2private(self, Integer3 coord):
        """global2private(g) -> Integer

        Transform a global coordinate to a private coordinate.

        Parameters
        ----------
        g : Integer3
            A global coordinate

        Returns
        -------
        Integer:
            A private coordinate

        """
        return self.thisptr.get().global2private(deref(coord.thisptr))

    def private2global(self, Integer coord):
        """private2global(coord) -> Integer3

        Transform a private coordinate to a global coordinate.

        """
        cdef Cpp_Integer3 g = self.thisptr.get().private2global(coord)
        return Integer3_from_Cpp_Integer3(address(g))

    def global2position(self, Integer3 g):
        """global2position(g) -> Real3

        Transform a global coordinate to a position.

        Parameters
        ----------
        g : Integer3
            A global coordinate

        Returns
        -------
        Real3:
            A position

        """
        cdef Cpp_Real3 pos = self.thisptr.get().global2position(deref(g.thisptr))
        return Real3_from_Cpp_Real3(address(pos))

    def position2global(self, Real3 pos):
        """position2global(pos) -> Integer3

        Transform a position to a global coordinate.

        Parameters
        ----------
        pos : Real3
            A position

        Returns
        -------
        Integer3:
            A global coordinate

        """
        cdef Cpp_Integer3 g = self.thisptr.get().position2global(deref(pos.thisptr))
        return Integer3_from_Cpp_Integer3(address(g))

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

    def as_base(self):
        """Return self as a base class. Only for developmental use."""
        retval = Space()
        del retval.thisptr
        retval.thisptr = new shared_ptr[Cpp_Space](
            <shared_ptr[Cpp_Space]>deref(self.thisptr))
        return retval

cdef LatticeWorld LatticeWorld_from_Cpp_LatticeWorld(
    shared_ptr[Cpp_LatticeWorld] w):
    r = LatticeWorld(Real3(1, 1, 1))
    r.thisptr.swap(w)
    return r

def create_lattice_world_cell_list_impl(
    edge_lengths, voxel_radius, matrix_sizes, rng):
    cdef shared_ptr[Cpp_LatticeWorld]* w = new shared_ptr[Cpp_LatticeWorld](
        create_lattice_world_cell_list_impl_alias(
            deref((<Real3>edge_lengths).thisptr), <Real>voxel_radius,
            deref((<Integer3>matrix_sizes).thisptr),
            deref((<GSLRandomNumberGenerator>rng).thisptr)))
    return LatticeWorld_from_Cpp_LatticeWorld(deref(w))

def create_lattice_world_vector_impl(edge_lengths, voxel_radius, rng):
    cdef shared_ptr[Cpp_LatticeWorld]* w = new shared_ptr[Cpp_LatticeWorld](
        create_lattice_world_vector_impl_alias(
            deref((<Real3>edge_lengths).thisptr), <Real>voxel_radius,
            deref((<GSLRandomNumberGenerator>rng).thisptr)))
    return LatticeWorld_from_Cpp_LatticeWorld(deref(w))

## LatticeSimulator
#  a python wrapper for Cpp_LatticeSimulator
cdef class LatticeSimulator:
    """ A class running the simulation with the lattice algorithm.

    LatticeSimulator(m, w)

    """

    def __init__(self, m, LatticeWorld w=None):
        """LatticeSimulator(m, w)
        LatticeSimulator(w)

        Constructor.

        Parameters
        ----------
        m : Model
            A model
        w : LatticeWorld
            A world

        """
        pass

    def __cinit__(self, m, LatticeWorld w=None):
        if w is None:
            self.thisptr = new Cpp_LatticeSimulator(
                deref((<LatticeWorld>m).thisptr))
        else:
            self.thisptr = new Cpp_LatticeSimulator(
                deref(Cpp_Model_from_Model(m)), deref(w.thisptr))

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
        cdef vector[pair[Cpp_ReactionRule, Cpp_ReactionInfo]] reactions = self.thisptr.last_reactions()
        cdef vector[pair[Cpp_ReactionRule, Cpp_ReactionInfo]].iterator it = reactions.begin()
        retval = []
        while it != reactions.end():
            retval.append((
                ReactionRule_from_Cpp_ReactionRule(
                    <Cpp_ReactionRule*>(address(deref(it).first))),
                ReactionInfo_from_Cpp_ReactionInfo(
                    <Cpp_ReactionInfo*>(address(deref(it).second)))))
            inc(it)
        return retval

    def set_alpha(self, Real alpha):
        """set_alpha(alpha)

        Set the value of alpha.

        Parameters
        ----------
        alpha : Real
            The value of alpha

        """
        self.thisptr.set_alpha(alpha)

    def get_alpha(self):
        """Return the value of alpha."""
        return self.thisptr.get_alpha()

    def calculate_alpha(self, ReactionRule rule):
        """calculate_alpha(rule) -> Real

        Return the recommended value of alpha

        Parameters
        ----------
        rule : ReactionRule
            A reaction rule.

        Returns
        -------
        Real:
            The recommneded value of alpha

        """
        return self.thisptr.calculate_alpha(deref(rule.thisptr))

    def model(self):
        """Return the model bound."""
        return Model_from_Cpp_Model(self.thisptr.model())

    def world(self):
        """Return the world bound."""
        return LatticeWorld_from_Cpp_LatticeWorld(self.thisptr.world())

    def run(self, Real duration, observers=None):
        """run(duration, observers)

        Run the simulation.

        Parameters
        ----------
        duration : Real
            A duration for running a simulation.
            A simulation is expected to be stopped at t() + duration.
        observers : list of Obeservers, optional
            observers

        """
        cdef vector[shared_ptr[Cpp_Observer]] tmp

        if observers is None:
            self.thisptr.run(duration)
        elif isinstance(observers, collections.Iterable):
            for obs in observers:
                tmp.push_back(deref((<Observer>(obs.as_base())).thisptr))
            self.thisptr.run(duration, tmp)
        else:
            self.thisptr.run(duration,
                deref((<Observer>(observers.as_base())).thisptr))

cdef LatticeSimulator LatticeSimulator_from_Cpp_LatticeSimulator(Cpp_LatticeSimulator* s):
    r = LatticeSimulator(
        Model_from_Cpp_Model(s.model()), LatticeWorld_from_Cpp_LatticeWorld(s.world()))
    del r.thisptr
    r.thisptr = s
    return r

## LatticeFactory
#  a python wrapper for Cpp_LatticeFactory
cdef class LatticeFactory:
    """ A factory class creating a LatticeWorld instance and a LatticeSimulator instance.

    LatticeFactory(Real voxel_radius=None, GSLRandomNumberGenerator rng=None)

    """

    def __init__(self, voxel_radius=None, GSLRandomNumberGenerator rng=None):
        """Constructor.

        Parameters
        ----------
        voxel_radius : Real, optional
            A radius of a voxel.
        rng : GSLRandomNumberGenerator, optional
            A random number generator.

        """
        pass

    def __cinit__(self, voxel_radius=None, GSLRandomNumberGenerator rng=None):
        if rng is not None:
            self.thisptr = new Cpp_LatticeFactory(<Real>voxel_radius, deref(rng.thisptr))
        elif voxel_radius is not None:
            self.thisptr = new Cpp_LatticeFactory(<Real>voxel_radius)
        else:
            self.thisptr = new Cpp_LatticeFactory()

    def __dealloc__(self):
        del self.thisptr

    def create_world(self, arg1=None):
        """create_world(arg1=None) -> LatticeWorld

        Return a LatticeWorld instance.

        Parameters
        ----------
        arg1 : Real3
            The lengths of edges of a LatticeWorld created

        or

        arg1 : str
            The path of a HDF5 file for LatticeWorld

        Returns
        -------
        LatticeWorld:
            The created world

        """
        if arg1 is None:
            return LatticeWorld_from_Cpp_LatticeWorld(
                shared_ptr[Cpp_LatticeWorld](self.thisptr.create_world()))
        elif isinstance(arg1, Real3):
            return LatticeWorld_from_Cpp_LatticeWorld(
                shared_ptr[Cpp_LatticeWorld](
                    self.thisptr.create_world(deref((<Real3>arg1).thisptr))))
        elif isinstance(arg1, str):
            return LatticeWorld_from_Cpp_LatticeWorld(
                shared_ptr[Cpp_LatticeWorld](self.thisptr.create_world(<string>(arg1))))
        else:
            return LatticeWorld_from_Cpp_LatticeWorld(
                shared_ptr[Cpp_LatticeWorld](self.thisptr.create_world(
                    deref(Cpp_Model_from_Model(arg1)))))

    def create_simulator(self, arg1, LatticeWorld arg2=None):
        """create_simulator(arg1, arg2) -> LatticeSimulator

        Return a LatticeSimulator instance.

        Parameters
        ----------
        arg1 : LatticeWorld
            A world

        or

        arg1 : Model
            A simulation model
        arg2 : LatticeWorld
            A world

        Returns
        -------
        LatticeSimulator:
            The created simulator

        """
        if arg2 is None:
            return LatticeSimulator_from_Cpp_LatticeSimulator(
                self.thisptr.create_simulator(deref((<LatticeWorld>arg1).thisptr)))
        else:
            return LatticeSimulator_from_Cpp_LatticeSimulator(
                self.thisptr.create_simulator(
                    deref(Cpp_Model_from_Model(arg1)), deref(arg2.thisptr)))
