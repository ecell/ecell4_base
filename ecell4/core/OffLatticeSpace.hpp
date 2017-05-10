#ifndef ECELL4_OFFLATTICE_SPACE_HPP
#define ECELL4_OFFLATTICE_SPACE_HPP

#include "VoxelSpaceBase.hpp"

namespace ecell4
{

class OffLatticeSpace : public VoxelSpaceBase
{
protected:

    typedef VoxelSpaceBase base_type;
    typedef std::vector<VoxelPool*> voxel_container;
    typedef std::vector<std::vector<coordinate_type> > adjoining_container;

public:
    typedef std::pair<coordinate_type, coordinate_type> coordinate_pair_type;
    typedef std::vector<Real3> position_container;
    typedef std::vector<coordinate_pair_type> coordinate_pair_list_type;

    OffLatticeSpace(const Real& voxel_radius);
    OffLatticeSpace(const Real& voxel_radius,
                    const position_container& positions,
                    const coordinate_pair_list_type& adjoining_pairs);
    virtual ~OffLatticeSpace();

    virtual std::pair<ParticleID, Voxel> get_voxel_at(const coordinate_type& coord) const;

    virtual const Particle particle_at(const coordinate_type& coord) const;

    virtual bool update_voxel(const ParticleID& pid, const Voxel& v);
    virtual bool remove_voxel(const ParticleID& pid);
    virtual bool remove_voxel(const coordinate_type& coord);

    virtual bool can_move(const coordinate_type& src,
                          const coordinate_type& dest) const;
    virtual bool move(const coordinate_type& src,
                      const coordinate_type& dest,
                      const std::size_t candidate=0);
    virtual std::pair<coordinate_type, bool> move_to_neighbor(
            VoxelPool* const& from,
            VoxelPool* const& loc,
            coordinate_id_pair_type& info,
            const Integer nrand);

    virtual VoxelPool* get_voxel_pool_at(const coordinate_type& coord) const;

    /*
     * Structure
     */
    virtual bool on_structure(const Voxel& v);

    virtual coordinate_type inner2coordinate(const coordinate_type inner) const;

    virtual Real3 coordinate2position(const coordinate_type& coord) const;
    virtual coordinate_type position2coordinate(const Real3& pos) const;

    virtual Integer num_neighbors(const coordinate_type& coord) const;
    virtual coordinate_type get_neighbor(const coordinate_type& coord,
                                         const Integer& nrand) const;
    virtual coordinate_type get_neighbor_boundary(const coordinate_type& coord,
                                                  const Integer& nrand) const;

    virtual Integer num_molecules(const Species& sp) const;

    virtual Real3 actual_lengths() const;

    virtual Integer size() const;
    virtual Integer3 shape() const;
    virtual Integer inner_size() const;

#ifdef WITH_HDF5
    virtual void save_hdf5(H5::Group* root) const;
    virtual void load_hdf5(const H5::Group& root);
#endif

protected:

    void reset(const position_container& positions,
               const coordinate_pair_list_type& adjoining_pairs);
    bool is_in_range(const coordinate_type& coord) const;
    VoxelPool* get_voxel_pool(const Voxel& v);
    coordinate_type get_coord(const ParticleID& pid) const;
    bool make_molecular_pool(const Species& sp,
                             Real radius,
                             Real D,
                             const std::string loc);
    Integer count_voxels(const boost::shared_ptr<VoxelPool>& vp) const;

protected:

    voxel_container voxels_;
    position_container positions_;
    adjoining_container adjoinings_;

    VoxelPool* vacant_;
};

} // ecell4

#endif /* ECELL4_OFFLATTICE_SPACE_HPP */
