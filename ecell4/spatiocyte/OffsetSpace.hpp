#ifndef ECELL4_SPATIOCYTE_OFFSET_SPACE_HPP
#define ECELL4_SPATIOCYTE_OFFSET_SPACE_HPP

#include <boost/shared_ptr.hpp>
#include <ecell4/core/Identifier.hpp>
#include <ecell4/core/Particle.hpp>
#include <ecell4/core/Real3.hpp>
#include <ecell4/core/Species.hpp>
#include <ecell4/core/ParticleVoxel.hpp>

namespace ecell4
{

namespace spatiocyte
{

template <typename T>
class OffsetSpace
{
protected:
    typedef typename T::coordinate_type coord_type;
    typedef std::pair<ParticleID, Particle> identified_particle;
    typedef std::pair<ParticleID, ParticleVoxel>    identified_voxel;

    boost::shared_ptr<T> space_;
    coord_type           offset_;

public:

    OffsetSpace(T* space, const coord_type& offset)
        : space_(space), offset_(offset)
    {}

    // OffsetSpace(boost::shared_ptr<T> space, const coord_type& offset)
    //     : space_(space), offset_(offset)
    // {}

    // operator OffsetSpace<const T>() const
    // {
    //     return OffsetSpace<const T>(space_, offset_);
    // }

    coord_type offset() const
    {
        return offset_;
    }

    void set_t(const Real& t)
    {
        space_->set_t(t);
    }

    bool is_in_range(const coord_type& coordinate) const
    {
        return coordinate >= offset_ && coordinate < offset_ + space_->size();
    }

    Real get_volume(const Species& species) const
    {
        return space_->get_volume(species);
    }

    std::size_t num_neighbors(const coord_type& coordinate) const
    {
        return space_->num_neighbors(coordinate - offset_);
    }

    coord_type get_neighbor(const coord_type& coordinate, std::size_t nrand) const
    {
        return space_->get_neighbor(coordinate - offset_, nrand) + offset_;
    }

    const Real3 coordinate2position(const coord_type& coordinate) const
    {
        return space_->coordinate2position(coordinate - offset_);
    }

    std::pair<ParticleID, ParticleVoxel> get_voxel(const ParticleID& pid) const
    {
        std::pair<ParticleID, ParticleVoxel> retval(space_->get_voxel(pid));
        retval.second.coordinate() -= offset_;
        return retval;
    }

    std::pair<ParticleID, ParticleVoxel> get_voxel_at(const coord_type& coordinate) const
    {
        return space_->get_voxel_at(coordinate - offset_);
    }

    boost::shared_ptr<const VoxelPool> get_voxel_pool_at(const coord_type& coordinate) const
    {
        return space_->get_voxel_pool_at(coordinate - offset_);
    }

    bool can_move(const coord_type& src, const coord_type& dest) const
    {
        return space_->can_move(src - offset_, dest - offset_);
    }

    //
    // Mutable functions
    //

    bool on_structure(ParticleVoxel voxel)
    {
        voxel.coordinate() -= offset_;
        return space_->on_structure(voxel);
    }

    bool update_voxel(const ParticleID& pid, ParticleVoxel voxel)
    {
        voxel.coordinate() -= offset_;
        return space_->update_voxel(pid, voxel);
    }

    bool remove_voxel(const ParticleID& pid)
    {
        return space_->remove_voxel(pid);
    }

    bool remove_particle(const ParticleID& pid)
    {
        return space_->remove_particle(pid);
    }

    bool remove_voxel(const coord_type& coordinate)
    {
        return space_->remove_voxel(coordinate - offset_);
    }

    bool move(const coord_type& src, const coord_type& dest, const std::size_t candidate)
    {
        return space_->move(src - offset_, dest - offset_, candidate);
    }

#define GET_MACRO( func, _1, name, ... ) name
#define WRAP( args... ) GET_MACRO( args, WRAP_GETTER_WITH_ARG, WRAP_GETTER )( args )
#define WRAP_GETTER( func ) func() const { return space_->func(); }
#define WRAP_GETTER_WITH_ARG( func, arg_type )\
    func(const arg_type& arg) const { return space_->func(arg); }

    const Real WRAP( t );

    const Real WRAP( volume );

    bool WRAP( has_species, Species );
    std::size_t WRAP( num_species );
    std::vector<Species> WRAP( list_species );

    const Real3& WRAP( edge_lengths );
    Real3 WRAP( actual_lengths );

    bool WRAP( has_particle, ParticleID );
    identified_particle WRAP( get_particle, ParticleID );
    std::size_t WRAP( num_particles                );
    std::size_t WRAP( num_particles,       Species );
    std::size_t WRAP( num_particles_exact, Species );
    std::vector<identified_particle> WRAP( list_particles                );
    std::vector<identified_particle> WRAP( list_particles,       Species );
    std::vector<identified_particle> WRAP( list_particles_exact, Species );

    std::size_t WRAP( num_molecules,       Species );
    std::size_t WRAP( num_molecules_exact, Species );

    Real WRAP( get_value,       Species );
    Real WRAP( get_value_exact, Species );

    Real WRAP( voxel_radius );
    Real WRAP( voxel_volume );
    Real WRAP( actual_volume );
    Real WRAP( unit_area );
    boost::shared_ptr<VoxelPool> WRAP( vacant );

    bool WRAP( has_voxel, ParticleID );
    std::size_t WRAP( num_voxels                );
    std::size_t WRAP( num_voxels,       Species );
    std::size_t WRAP( num_voxels_exact, Species );
    std::vector<identified_voxel> WRAP( list_voxels                );
    std::vector<identified_voxel> WRAP( list_voxels,       Species );
    std::vector<identified_voxel> WRAP( list_voxels_exact, Species );

    bool WRAP( has_molecule_pool, Species );
    boost::shared_ptr<const VoxelPool>    WRAP( find_voxel_pool,    Species );
    boost::shared_ptr<const MoleculePool> WRAP( find_molecule_pool, Species );

    coord_type WRAP( position2coordinate, Real3 );
    const Integer WRAP( size );
    const Integer3 WRAP( shape );

#undef GET_MACRO
#undef WRAP
#undef WRAP_GETTER
#undef WRAP_GETTER_WITH_ARG

    boost::shared_ptr<VoxelPool> find_voxel_pool(const Species& species)
    {
        return space_->find_voxel_pool(species);
    }

    boost::shared_ptr<MoleculePool> find_molecule_pool(const Species& species)
    {
        return space_->find_molecule_pool(species);
    }

    bool
    make_structure_type(const Species& sp, Shape::dimension_kind dimension, const std::string loc)
    {
        return space_->make_structure_type(sp, dimension, loc);
    }

    bool
    make_interface_type(const Species& sp, Shape::dimension_kind dimension, const std::string loc)
    {
        return space_->make_interface_type(sp, dimension, loc);
    }

#ifdef WITH_HDF5
    void save_hdf5(H5::Group* root) const
    {
        space_->save_hdf5(root);
    }

    void load_hdf5(const H5::Group& root)
    {
        space_->load_hdf5(root);
    }
#endif

};

} // ns spatiocyte

} // ns ecell4

#endif // ECELL4_SPATIOCYTE_OFFSET_SPACE_HPP
