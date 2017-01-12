#ifndef __ECELL4_OFFLATTICE_SPACE_HPP
#define __ECELL4_OFFLATTICE_SPACE_HPP

#include "VoxelSpaceBase.hpp"

namespace ecell4
{

class OffLatticeSpace : public VoxelSpaceBase
{
protected:

    typedef VoxelSpaceBase base_type;
    typedef std::vector<VoxelPool*> voxel_container;
    typedef std::vector<std::vector<voxel_container::iterator> > adjoining_container;

public:
    typedef std::vector<Real3> position_container;

    OffLatticeSpace(const Real& voxel_radius, const position_container& positions);
    virtual ~OffLatticeSpace();

    virtual std::pair<ParticleID, Voxel> get_voxel_at(const coordinate_type& coord) const;

    virtual const Particle particle_at(const coordinate_type& coord) const;

    virtual bool update_voxel(const ParticleID& pid, const Voxel& v);
    virtual bool remove_voxel(const ParticleID& pid);
    virtual bool remove_voxel(const coordinate_type& coord);

    virtual bool move(
        const coordinate_type& src, const coordinate_type& dest,
        const std::size_t candidate=0);
    virtual std::pair<coordinate_type, bool> move_to_neighbor(
        VoxelPool* const& from, VoxelPool* const& loc,
        coordinate_id_pair_type& info, const Integer nrand);

    virtual VoxelPool* get_voxel_pool_at(const coordinate_type& coord) const;

    /*
     * Structure
     */
    virtual bool on_structure(const Voxel& v);

    virtual coordinate_type inner2coordinate(const coordinate_type inner) const;

    virtual Real3 coordinate2position(const coordinate_type& coord) const;
    virtual coordinate_type position2coordinate(const Real3& pos) const;

    virtual coordinate_type get_neighbor(
        const coordinate_type& coord, const Integer& nrand) const;
    virtual coordinate_type get_neighbor_boundary(
        const coordinate_type& coord, const Integer& nrand) const;

    virtual Integer num_molecules(const Species& sp) const;

    virtual Integer size() const;
    virtual Integer3 shape() const;
    virtual Integer inner_size() const;

protected:

    void reset(const position_container& positions);
    bool is_in_range(const coordinate_type& coord) const;
    VoxelPool* get_voxel_pool(const Voxel& v);
    coordinate_type get_coord(const ParticleID& pid) const;
    bool make_molecular_pool(const Species& sp, Real radius, Real D, const std::string loc);
    Integer count_voxels(const boost::shared_ptr<VoxelPool>& vp) const;

protected:

    voxel_container voxels_;
    position_container positions_;
    adjoining_container adjoinings_;

    VoxelPool* vacant_;
};

} // ecell4

#endif /* __ECELL4_OFFLATTICE_SPACE_HPP */
