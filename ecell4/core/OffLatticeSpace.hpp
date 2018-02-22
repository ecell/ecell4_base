#ifndef ECELL4_OFFLATTICE_SPACE_HPP
#define ECELL4_OFFLATTICE_SPACE_HPP

#include <boost/optional.hpp>
#include "VoxelSpaceBase.hpp"

namespace ecell4
{

class OffLatticeSpace : public VoxelSpaceBase
{
protected:

    typedef VoxelSpaceBase base_type;
    typedef std::vector<boost::shared_ptr<VoxelPool> > voxel_container;
    typedef std::vector<std::vector<coordinate_type> > adjoining_container;

public:
    typedef std::pair<coordinate_type, coordinate_type> coordinate_pair_type;
    typedef std::vector<Real3> position_container;
    typedef std::vector<coordinate_pair_type> coordinate_pair_list_type;

    /*
     * Constructor and Destructor
     */
    OffLatticeSpace(const Real& voxel_radius);
    OffLatticeSpace(const Real& voxel_radius,
                    const position_container& positions,
                    const coordinate_pair_list_type& adjoining_pairs);
    ~OffLatticeSpace();

    /*
     * Space Traits
     */
#ifdef WITH_HDF5
    void save_hdf5(H5::Group* root) const
    {
        throw NotSupported("OffLatticeSpace::save_hdf5 is not supported.");
    }

    void load_hdf5(const H5::Group& root)
    {
        throw NotSupported("OffLatticeSpace::load_hdf5 is not supported.");
    }
#endif

    /*
     * ParticleSpace Traits
     */
    Real3 actual_lengths() const
    {
        throw NotSupported("OffLatticeSpace::actual_lengths() is not supported.");
    }

    const Particle particle_at(const coordinate_type& coord) const;

    /*
     * VoxelSpace Traits
     */
    std::pair<ParticleID, Voxel> get_voxel_at(const coordinate_type& coord) const;
    boost::shared_ptr<VoxelPool> get_voxel_pool_at(const coordinate_type& coord) const
    {
        return voxels_.at(coord);
    }

    /*
     * Coordinate Transformation
     */
    coordinate_type inner2coordinate(const coordinate_type inner) const
    {
        throw NotSupported("OffLatticeSpace::inner2coordinate() is not supported.");
    }

    Real3 coordinate2position(const coordinate_type& coord) const
    {
        return positions_.at(coord);
    }

    coordinate_type position2coordinate(const Real3& pos) const;

    /*
     * Neighbor
     */
    Integer num_neighbors(const coordinate_type& coord) const
    {
        return adjoinings_.at(coord).size();
    }

    coordinate_type
    get_neighbor(const coordinate_type& coord, const Integer& nrand) const
    {
        return adjoinings_.at(coord).at(nrand);
    }

    /*
     * Voxel Manipulation
     */
    bool update_voxel(const ParticleID& pid, const Voxel& v);
    bool remove_voxel(const ParticleID& pid);
    bool remove_voxel(const coordinate_type& coord);

    bool can_move(const coordinate_type& src,
                  const coordinate_type& dest) const;

    bool move(const coordinate_type& src,
              const coordinate_type& dest,
              const std::size_t candidate=0);

    Integer size() const
    {
        return voxels_.size();
    }

    Integer3 shape() const
    {
        throw NotSupported("OffLatticeSpace::shape() is not supported.");
    }

    Integer inner_size() const
    {
        return size();
    }

protected:

    boost::shared_ptr<VoxelPool> get_voxel_pool(const Voxel& v);
    Integer count_voxels(const boost::shared_ptr<VoxelPool>& vp) const
    {
        return static_cast<Integer>(std::count(voxels_.begin(), voxels_.end(), vp));
    }

    bool is_in_range(const coordinate_type& coord) const
    {
        return 0 <= coord && coord < voxels_.size();
    }

    void reset(const position_container& positions,
               const coordinate_pair_list_type& adjoining_pairs);
    boost::optional<coordinate_type> get_coord(const ParticleID& pid) const;
    bool make_molecular_pool(const Species& sp,
                             Real radius,
                             Real D,
                             const std::string loc);

protected:

    voxel_container voxels_;
    position_container positions_;
    adjoining_container adjoinings_;

};

} // ecell4

#endif /* ECELL4_OFFLATTICE_SPACE_HPP */
