from libcpp.string cimport string
from libcpp cimport bool
from libcpp.vector cimport vector

from ecell4.types cimport *
from ecell4.shared_ptr cimport shared_ptr
from ecell4.core cimport *


## Cpp_SpatiocyteWorld
#  ecell4::spatiocyte::SpatiocyteWorld
cdef extern from "ecell4/spatiocyte/SpatiocyteWorld.hpp" namespace "ecell4::spatiocyte":
    cdef cppclass Cpp_SpatiocyteWorld "ecell4::spatiocyte::SpatiocyteWorld":
        Cpp_SpatiocyteWorld(Cpp_Position3& Position3, Real& voxel_radius) except +
        void set_t(Real t)
        Real t()
        Real voxel_radius()
        # boost::array<Integer, 3> lattice_size()
        Integer num_voxels()
        Cpp_Voxel get_voxel(Integer& coord)
        vector[unsigned int] coordinates(Cpp_Species& sp)
        Cpp_Position3 edge_lengths()
        Integer num_particles()
        Integer num_particles(Cpp_Species& sp)
        vector[pair[Cpp_ParticleID, Cpp_Particle]] list_particles(Cpp_Species& sp)
        Real volume()
        bool has_species(Cpp_Species& sp)
        Integer num_molecules(Cpp_Species& sp)
        void add_molecules(Cpp_Species& sp, Integer num)
        Real dt()
        void save(string filename)

## SpatiocyteWorld
#  a python wrapper for Cpp_SpatiocyteWorld
cdef class SpatiocyteWorld:
    cdef shared_ptr[Cpp_SpatiocyteWorld]* thisptr

## Cpp_SpatiocyteSimulator
#  ecell4::spatiocyte::SpatiocyteSimulator
cdef extern from "ecell4/spatiocyte/SpatiocyteSimulator.hpp" namespace "ecell4::spatiocyte":
    cdef cppclass Cpp_SpatiocyteSimulator "ecell4::spatiocyte::SpatiocyteSimulator":
        Cpp_SpatiocyteSimulator(
            shared_ptr[Cpp_NetworkModel], shared_ptr[Cpp_SpatiocyteWorld]) except +
        Integer num_steps()
        void step()
        bool step(Real& upto)
        Real t()
        Real dt()
        void initialize()

## SpatiocyteSimulator
#  a python wrapper for Cpp_SpatiocyteSimulator
cdef class SpatiocyteSimulator:
    cdef Cpp_SpatiocyteSimulator* thisptr
