import collections
from cython.operator cimport dereference as deref, preincrement as inc
from cython cimport address
from libcpp.string cimport string
from libcpp.vector cimport vector

from ecell4.types cimport *
from ecell4.shared_ptr cimport shared_ptr
from ecell4.core cimport *

## LatticeWorld
#  a python wrapper for Cpp_LatticeWorld
cdef class LatticeWorld:
    """A class containing the properties of the simulation world

    LatticeWorld(edge_lengths=None, voxel_radius=None, GSLRandomNumberGenerator rng=None)

    """

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
        """Set the value of the time of the world
        Args:
            t (Real): the time of the world
        """
        self.thisptr.get().set_t(t)

    def t(self):
        """Return the time of the world
        Returns:
            Real: the time of the world
        """
        return self.thisptr.get().t()

    def volume(self):
        """Return the volume of the world
        Returns:
            Real: the volume of the world
        """
        return self.thisptr.get().volume()

    def voxel_volume(self):
        """Return the volume of a voxel
        Returns:
            Real: the volume of a voxel
        """
        return self.thisptr.get().voxel_volume()

    def get_volume(self):
        """Return the actual volume of the world
        Returns:
            Real: the actual volume of the world
        """
        return self.thisptr.get().get_volume()

    def actual_lengths(self):
        """Return the actual edge lengths of the world
        Returns:
            Real3: the actual edge lengths of the world
        """
        cdef Cpp_Real3 lengths = self.thisptr.get().actual_lengths()
        return Real3_from_Cpp_Real3(address(lengths))

    def new_particle(self, arg1, Real3 arg2=None):
        """Create a new particle
        Args:
            arg1 (Particle): a particle container
        or
        Args:
            arg1 (Species): a species of a particle
            arg2 (Real3): a coordinate to place a particle
        Returns:
            tuple: a pair of ParticleID and Particle of a new particle
        """
        cdef pair[pair[Cpp_ParticleID, Cpp_Particle], bool] retval

        if arg2 is None:
            retval = self.thisptr.get().new_particle(deref((<Particle> arg1).thisptr))
        else:
            retval = self.thisptr.get().new_particle(deref((<Species> arg1).thisptr), deref(arg2.thisptr))
        return ((ParticleID_from_Cpp_ParticleID(address(retval.first.first)), Particle_from_Cpp_Particle(address(retval.first.second))), retval.second)

    def get_particle(self, ParticleID pid):
        """Return the particle associated a given ParticleID
        Args:
            pid (ParticleID): a id of the particle you want
        Returns:
            tuple: a pair of ParticleID and Particle
        """
        cdef pair[Cpp_ParticleID, Cpp_Particle] \
            pid_particle_pair = self.thisptr.get().get_particle(deref(pid.thisptr))
        return (ParticleID_from_Cpp_ParticleID(address(pid_particle_pair.first)),
                Particle_from_Cpp_Particle(address(pid_particle_pair.second)))

    def get_voxel(self, arg):
        """Return the voxel having a particle associated with a given ParticleID or coordinate
        Args:
            arg (ParticleID or Integer): an id or coordiante of the particle in the voxel you want
        Returns:
            tuple: a pair of ParticleID and Voxel
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
        """Check if the given species would be on the proper structure at the coordinate
        Args:
            sp (Species): a species scheduled to be placed
            coord (Integer): a coordinate to be occupied
        Returns:
            bool: if it is on the proper structure, or not
        """
        return self.thisptr.get().on_structure(deref(sp.thisptr), coord)

    def remove_particle(self, ParticleID pid):
        """Remove the particle associated with a given ParticleID
        Args:
            pid (ParticleID): a id of particle to remove
        """
        self.thisptr.get().remove_particle(deref(pid.thisptr))

    def remove_voxel(self, ParticleID pid):
        """Remove the particle associated with a given ParticleID
        Args:
            pid (ParticleID): a id of particle to remove
        """
        self.thisptr.get().remove_voxel(deref(pid.thisptr))

    def edge_lengths(self):
        """Return the edge lengths of the world
        Returns:
            Real3: the edge lengths of the world
        """
        cdef Cpp_Real3 lengths = self.thisptr.get().edge_lengths()
        return Real3_from_Cpp_Real3(address(lengths))

    def num_particles(self, Species sp = None):
        """Return the number of particles
        Args:
            sp (Species, optional): the species of particles to count
        Returns:
            Integer: the number of particles (of the given species)
        """
        if sp is None:
            return self.thisptr.get().num_particles()
        else:
            return self.thisptr.get().num_particles(deref(sp.thisptr))

    def num_particles_exact(self, Species sp):
        """Return the number of particles of a given species
        Args:
            sp (Species): the species of particles to count
        Returns:
            Integer: the number of particles of a given species
        """
        return self.thisptr.get().num_particles_exact(deref(sp.thisptr))

    def num_voxels(self, Species sp = None):
        """Return the number of voxels
        Args:
            sp (Species, optional): the species of particles to count
        Returns:
            Integer: the number of voxels (of the given species)
        """
        if sp is None:
            return self.thisptr.get().num_voxels()
        else:
            return self.thisptr.get().num_voxels(deref(sp.thisptr))

    def num_voxels_exact(self, Species sp):
        """Return the number of voxels of a given species
        Args:
            sp (Species): the species of particles to count
        Returns:
            Integer: the number of voxels of a given species
        """
        return self.thisptr.get().num_voxels_exact(deref(sp.thisptr))

    def list_particles(self, Species sp = None):
        """Return the list of particles
        Args:
            sp (Species, optional): the species of particles to list up
        Returns:
            list: the list of particles (of the given species)
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
        """Return the list of particles of a given species
        Args:
            sp (Species): the species of particles to list up
        Returns:
            list: the list of particles of a given species
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
        """Return the neighbor voxel of a given voxel
        Args:
            coord (Integer): a coordinate of a voxel
            nrand (Integer): a key in the range from 0 to 11 to assign a neighbor voxel
        Returns:
            Integer: the coordinate of the neighbor voxel
        """
        return self.thisptr.get().get_neighbor(coord, nrand)

    def get_neighbor_private(self, coord, nrand):
        """Return the neighbor voxel of a given voxel
        Args:
            coord (Integer): a private coordinate of a voxel
            nrand (Integer): a key in the range from 0 to 11 to assign a neighbor voxel
        Returns:
            Integer: the private coordinate of the neighbor voxel
        """
        return self.thisptr.get().get_neighbor_private(coord, nrand)

    def has_particle(self, ParticleID pid):
        """Check if a particle associated with a given particle id exists
        Args:
            pid (ParticleID): a particle id to check
        Returns:
            bool: if a particle exists, this is true. Otherwise false
        """
        return self.thisptr.get().has_particle(deref(pid.thisptr))

    def update_particle(self, ParticleID pid, Particle p):
        """Update a particle
        Args:
            pid (ParticleID): a particle id of the particle to update
            p (Particle): the information to update a particle
        """
        return self.thisptr.get().update_particle(deref(pid.thisptr), deref(p.thisptr))

    def num_molecules(self, Species sp = None):
        """Return the number of molecules
        Args:
            sp (Species, optional): a species whose molecules you count
        Returns:
            Integer: the number of molecules (of a given species)
        """
        if sp is None:
            return self.thisptr.get().num_molecules()
        else:
            return self.thisptr.get().num_molecules(deref(sp.thisptr))

    def num_molecules_exact(self, Species sp):
        """Return the number of molecules of a given species
        Args:
            sp (Species): a species whose molecules you count
        Returns:
            Integer: the number of molecules of a given species
        """
        return self.thisptr.get().num_molecules_exact(deref(sp.thisptr))

    def add_molecules(self, Species sp, Integer num, shape=None):
        """Add some molecules
        Args:
            sp (Species): a species of molecules to add
            num (Integer): the number of molecules to add
            shape (Shape, optional): a shape to add molecules on
        """
        if shape is None:
            self.thisptr.get().add_molecules(deref(sp.thisptr), num)
        else:
            self.thisptr.get().add_molecules(
                deref(sp.thisptr), num, deref((<Shape>(shape.as_base())).thisptr))

    def remove_molecules(self, Species sp, Integer num):
        """Remove the molecules
        Args:
            sp (Species): a species whose molecules to remove
        """
        self.thisptr.get().remove_molecules(deref(sp.thisptr), num)

    def save(self, filename):
        """Save the world to a file
        Args:
            filename (str): a filename to save to
        """
        self.thisptr.get().save(tostring(filename))

    def load(self, filename):
        """Load the world from a file
        Args:
            filename (str): a filename to load from
        """
        self.thisptr.get().load(tostring(filename))

    def new_voxel(self, arg1, arg2=None):
        """Create a particle
        Args:
            arg1 (Voxel): the information to create
        or
        Args:
            arg1 (Species): the Species of particles to create
            arg2 (Integer): the number of particles(voxels)
        Returns:
            tuple: a pair of ParticleID and Voxel
        """
        cdef pair[pair[Cpp_ParticleID, Cpp_Voxel], bool] retval

        if arg2 is None:
            retval = self.thisptr.get().new_voxel(deref((<Voxel> arg1).thisptr))
        else:
            retval = self.thisptr.get().new_voxel(deref((<Species> arg1).thisptr), <Integer> arg2)
        return ((ParticleID_from_Cpp_ParticleID(address(retval.first.first)), Voxel_from_Cpp_Voxel(address(retval.first.second))), retval.second)

    def update_voxel(self, ParticleID pid, Voxel v):
        """Update a particle
        Args:
            pid (ParticleID): a particle id of the particle to update
            v (Voxel): the information to update
        Returns:
            bool: whether to succeed to update the particle
        """
        return self.thisptr.get().update_voxel(deref(pid.thisptr), deref(v.thisptr))

    def list_voxels(self, Species sp = None):
        """Returns the list of voxels
        Args:
            sp (Species, optional): a species of particles to list up
        Returns:
            list: the list of the pair of ParticleID and Voxel
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
        """Check if a particle exists
        Args:
            pid (ParticleID): a particle id of the particle to check
        Returns:
            bool: whether a particle associated with a given particle id exists
        """
        return self.thisptr.get().has_voxel(deref(pid.thisptr))

    def voxel_radius(self):
        """Return the voxel radius
        Returns:
            Real: the voxel radius
        """
        return self.thisptr.get().voxel_radius()

    def col_size(self):
        """Return the size of the column of the world
        Returns:
            Integer: the size of the column of the world
        """
        return self.thisptr.get().col_size()

    def row_size(self):
        """Return the size of row of the world
        Returns:
            Integer: the size of the row of the world
        """
        return self.thisptr.get().row_size()

    def layer_size(self):
        """Return the size of layer of the world
        Returns:
            Integer: the size of layer of the world
        """
        return self.thisptr.get().layer_size()

    def size(self):
        """Return the size of voxels
        Returns:
            Integer: the size of voxels
        """
        return self.thisptr.get().size()

    def shape(self):
        """Return the triplet of sizes of column, row and layer
        Returns:
            Integer3: the triplet of sizes of column, row and layer
        """
        cdef Cpp_Integer3 sizes = self.thisptr.get().shape()
        return Integer3_from_Cpp_Integer3(address(sizes))

    def bind_to(self, m):
        """Bind a model to the world
        Args:
            m (Model): a model to bind
        """
        self.thisptr.get().bind_to(deref(Cpp_Model_from_Model(m)))

    def private2position(self, Integer coord):
        """Transform private coordinate to position
        """
        cdef Cpp_Real3 pos = self.thisptr.get().private2position(coord)
        return Real3_from_Cpp_Real3(address(pos))

    def coordinate2position(self, Integer coord):
        cdef Cpp_Real3 pos = self.thisptr.get().coordinate2position(coord)
        return Real3_from_Cpp_Real3(address(pos))

    def position2coordinate(self, Real3 pos):
        return self.thisptr.get().position2coordinate(
            deref(pos.thisptr))

    def private2coord(self, Integer coord):
        return self.thisptr.get().private2coord(coord)

    def coord2private(self, Integer coord):
        return self.thisptr.get().coord2private(coord)

    def global2coord(self, Integer3 coord):
        return self.thisptr.get().global2coord(deref(coord.thisptr))

    def coord2global(self, Integer coord):
        cdef Cpp_Integer3 g = self.thisptr.get().coord2global(coord)
        return Integer3_from_Cpp_Integer3(address(g))

    def global2private(self, Integer3 coord):
        return self.thisptr.get().global2private(deref(coord.thisptr))

    def private2global(self, Integer coord):
        cdef Cpp_Integer3 g = self.thisptr.get().private2global(coord)
        return Integer3_from_Cpp_Integer3(address(g))

    def global2position(self, Integer3 g):
        cdef Cpp_Real3 pos = self.thisptr.get().global2position(deref(g.thisptr))
        return Real3_from_Cpp_Real3(address(pos))

    def position2global(self, Real3 pos):
        cdef Cpp_Integer3 g = self.thisptr.get().position2global(deref(pos.thisptr))
        return Integer3_from_Cpp_Integer3(address(g))

    def add_structure(self, Species sp, shape):
        return self.thisptr.get().add_structure(
            deref(sp.thisptr), deref((<Shape>(shape.as_base())).thisptr))

    def rng(self):
        return GSLRandomNumberGenerator_from_Cpp_RandomNumberGenerator(
            self.thisptr.get().rng())

    def as_base(self):
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
    """ A class running the simulation according to the lattice model.

    LatticeSimulator(m, LatticeWorld w=None)

    """

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
        """Return the number of steps
        Returns:
            Integer: the number of steps
        """
        return self.thisptr.num_steps()

    def step(self, upto = None):
        """Step the simulation
        Args:
            upto (Real): the time which to step the simulation up to
        """
        if upto is None:
            self.thisptr.step()
        else:
            return self.thisptr.step(upto)

    def t(self):
        """Return the time
        Returns:
            Real: the time
        """
        return self.thisptr.t()

    def dt(self):
        """Return the interval of times
        Returns:
            Real: the interval of times
        """
        return self.thisptr.dt()

    def next_time(self):
        """Return the time of the next step
        Returns:
            Real: the time of the next step
        """
        return self.thisptr.next_time()

    def set_dt(self, Real dt):
        """Set the interval of times
        Args:
            dt (Real): the interval of times
        """
        self.thisptr.set_dt(dt)

    def initialize(self):
        """Initialize the simulator
        """
        self.thisptr.initialize()

    def last_reactions(self):
        """Return reactions occuring at the last step
        Returns:
            list: the list of reactions
        """
        cdef vector[Cpp_ReactionRule] reactions = self.thisptr.last_reactions()
        cdef vector[Cpp_ReactionRule].iterator it = reactions.begin()
        retval = []
        while it != reactions.end():
            retval.append(ReactionRule_from_Cpp_ReactionRule(
                <Cpp_ReactionRule*>(address(deref(it)))))
            inc(it)
        return retval

    def set_alpha(self, Real alpha):
        """Set the value of alpha
        Args:
            alpha (Real): the value of alpha
        """
        self.thisptr.set_alpha(alpha)

    def get_alpha(self):
        """Return the value of alpha
        Returns:
            Real: the value of alpha
        """
        return self.thisptr.get_alpha()

    def model(self):
        """Return the reaction model
        Returns:
            Model: the reaction model
        """
        return Model_from_Cpp_Model(self.thisptr.model())

    def world(self):
        """Return the world
        Returns:
            LatticeWorld: the world
        """
        return LatticeWorld_from_Cpp_LatticeWorld(self.thisptr.world())

    def run(self, Real duration, observers=None):
        """Run the simulation
        Args:
            duration (Real): the time in second
            observers (list of Obeservers, optional): the observers
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

    LatticeFactory(voxel_radius=None, GSLRandomNumberGenerator rng=None)

    """

    def __cinit__(self, voxel_radius=None, GSLRandomNumberGenerator rng=None):
        if rng is not None:
            self.thisptr = new Cpp_LatticeFactory(<Real>voxel_radius, deref(rng.thisptr))
        elif voxel_radius is not None:
            self.thisptr = new Cpp_LatticeFactory(<Real>voxel_radius)
        else:
            self.thisptr = new Cpp_LatticeFactory()

    def __dealloc__(self):
        del self.thisptr

    def create_world(self, arg1):
        """Return a LatticeWorld instance.
        Args:
            arg1 (Real3): The lengths of edges of a LatticeWorld created
        or
        Args:
            arg1 (str): The path of a file of LatticeWorld
        Returns:
            LatticeWorld: the created world
        """
        if isinstance(arg1, Real3):
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
        """Return a LatticeSimulator instance.
        Args:
            arg1 (LatticeWorld): a world in which the simulation runs
        or
        Args:
            arg1 (Model): a simulation model
            arg2 (LatticeWorld): a world in which the simulation runs
        Returns:
            LatticeSimulator: the created simulator
        """
        if arg2 is None:
            return LatticeSimulator_from_Cpp_LatticeSimulator(
                self.thisptr.create_simulator(deref((<LatticeWorld>arg1).thisptr)))
        else:
            return LatticeSimulator_from_Cpp_LatticeSimulator(
                self.thisptr.create_simulator(
                    deref(Cpp_Model_from_Model(arg1)), deref(arg2.thisptr)))
