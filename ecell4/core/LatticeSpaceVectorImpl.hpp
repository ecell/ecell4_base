#ifndef ECELL4_LATTICE_SPACE_VECTOR_IMPL_HPP
#define ECELL4_LATTICE_SPACE_VECTOR_IMPL_HPP

#include "HCPLatticeSpace.hpp"

namespace ecell4 {

class LatticeSpaceVectorImpl
    : public HCPLatticeSpace
{
public:

    typedef HCPLatticeSpace base_type;
    typedef std::vector<boost::shared_ptr<VoxelPool> > voxel_container;

public:

    LatticeSpaceVectorImpl(const Real3& edge_lengths,
                           const Real& voxel_radius,
                           const bool is_periodic = true);
    ~LatticeSpaceVectorImpl();

    /*
     * Space APIs
     *
     * using ParticleID, Species and Posision3
     */

    Integer num_species() const;

    bool remove_voxel(const ParticleID& pid);
    bool remove_voxel(const coordinate_type& coord);

    bool update_structure(const Particle& p);

    /*
     * for Simulator
     *
     * using Species and coordinate_type
     */
    std::vector<std::pair<ParticleID, ParticleVoxel> > list_voxels() const;
    std::vector<std::pair<ParticleID, ParticleVoxel> > list_voxels(const Species& sp) const;
    std::vector<std::pair<ParticleID, ParticleVoxel> > list_voxels_exact(const Species& sp) const;

    std::pair<ParticleID, ParticleVoxel> get_voxel_at(const coordinate_type& coord) const;

    bool update_voxel(const ParticleID& pid, ParticleVoxel v);
    bool add_voxel(const Species& species, const ParticleID& pid, const coordinate_type& coord);

    bool add_voxels(const Species& species,
                    std::vector<std::pair<ParticleID, coordinate_type> > voxels);

    const Species& find_species(std::string name) const;
    std::vector<coordinate_type> list_coords(const Species& sp) const;
    std::vector<coordinate_type> list_coords_exact(const Species& sp) const;

    boost::shared_ptr<VoxelPool> get_voxel_pool_at(const coordinate_type& coord) const
    {
        return voxels_.at(coord);
    }

    bool move(const coordinate_type& src,
              const coordinate_type& dest,
              const std::size_t candidate=0);
    bool can_move(const coordinate_type& src, const coordinate_type& dest) const;

    coordinate_type
    get_neighbor(const coordinate_type& coord, const Integer& nrand) const
    {
        coordinate_type const dest = get_neighbor_(coord, nrand);

        if (voxels_.at(dest) != periodic_)
        {
            return dest;
        }
        else
        {
            return periodic_transpose(dest);
        }
    }

    bool is_periodic() const
    {
        return is_periodic_;
    }

#ifdef WITH_HDF5
    /*
     * HDF5 Save
     */
    void save_hdf5(H5::Group* root) const
    {
        save_lattice_space(*this, root, "LatticeSpaceVectorImpl");
    }

    void load_hdf5(const H5::Group& root)
    {
        load_lattice_space(root, this);
    }
#endif

    void reset(const Real3& edge_lengths, const Real& voxel_radius, const bool is_periodic)
    {
        base_type::reset(edge_lengths, voxel_radius, is_periodic);

        is_periodic_ = is_periodic;
        initialize_voxels(is_periodic_);
    }

protected:

    coordinate_type apply_boundary_(const coordinate_type& coord) const
    {
        return periodic_transpose(coord);
    }

    void initialize_voxels(const bool is_periodic);

    std::pair<coordinate_type, bool>
    move_(coordinate_type from,
          coordinate_type to,
          const std::size_t candidate=0);

    std::pair<coordinate_type, bool>
    move_(coordinate_id_pair_type& info,
          coordinate_type to);

    coordinate_type get_coord(const ParticleID& pid) const;

protected:

    bool is_periodic_;

    voxel_container voxels_;

    boost::shared_ptr<VoxelPool> border_;
    boost::shared_ptr<VoxelPool> periodic_;
};

} // ecell4

#endif
