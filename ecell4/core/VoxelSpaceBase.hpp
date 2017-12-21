#ifndef ECELL4_VOXELSPACEBASE_HPP
#define ECELL4_VOXELSPACEBASE_HPP

#include "LatticeSpace.hpp"

namespace ecell4
{

static inline std::string get_location_serial(const VoxelPool* vp)
{
    if (vp == NULL || vp->location() == NULL || vp->location()->is_vacant()) {
        return "";
    }

    return vp->location()->species().serial();
}

static inline std::string get_location_serial(const boost::shared_ptr<VoxelPool>& vp)
{
    if (vp == NULL || vp->location() == NULL || vp->location()->is_vacant()) {
        return "";
    }

    return vp->location()->species().serial();
}


class VoxelSpaceBase : public LatticeSpace
{
protected:

    typedef LatticeSpace base_type;
    typedef utils::get_mapper_mf<Species, boost::shared_ptr<VoxelPool> >::type
        voxel_pool_map_type;
    typedef utils::get_mapper_mf<Species, boost::shared_ptr<MoleculePool> >::type
        molecule_pool_map_type;

public:

    VoxelSpaceBase(const Real& voxel_radius);
    virtual ~VoxelSpaceBase();

    const Real t() const
    {
        return t_;
    }

    void set_t(const Real& t)
    {
        if (t < 0.0)
        {
            throw std::invalid_argument("the time must be positive.");
        }
        t_ = t;
    }

    Real voxel_radius() const
    {
        return voxel_radius_;
    }

    Real voxel_volume() const
    {
        return calculate_voxel_volume(voxel_radius_);
    }

    Real get_volume(const Species& sp) const
    {
        return voxel_volume() * num_voxels_exact(sp);
        // return inner_size() * voxel_volume();
    }

    Real actual_volume() const
    {
        return inner_size() * voxel_volume();
    }

    Real unit_area() const
    {
        const Real r(voxel_radius_);
        return 2.0 * sqrt(3.0) * r * r;
    }


    std::vector<Species> list_species() const;

    Integer num_voxels_exact(const Species& sp) const;
    Integer num_voxels(const Species& sp) const;
    Integer num_voxels() const;
    bool has_voxel(const ParticleID& pid) const;

    std::vector<std::pair<ParticleID, Voxel> > list_voxels() const;
    std::vector<std::pair<ParticleID, Voxel> > list_voxels(const Species& sp) const;
    std::vector<std::pair<ParticleID, Voxel> > list_voxels_exact(const Species& sp) const;

    std::pair<ParticleID, Voxel> get_voxel(const ParticleID& pid) const;
    VoxelPool* find_voxel_pool(const Species& sp);
    const VoxelPool* find_voxel_pool(const Species& sp) const;
    bool has_molecule_pool(const Species& sp) const;
    MoleculePool* find_molecule_pool(const Species& sp);
    const MoleculePool* find_molecule_pool(const Species& sp) const;

    /*
     * Virtual functions
     */

    virtual Real3 actual_lengths() const = 0;

    virtual std::pair<ParticleID, Voxel> get_voxel_at(const coordinate_type& coord) const = 0;
    virtual VoxelPool* get_voxel_pool_at(const coordinate_type& coord) const = 0;
    virtual const Particle particle_at(const coordinate_type& coord) const = 0;

    virtual bool update_voxel(const ParticleID& pid, const Voxel& v) = 0;
    virtual bool remove_voxel(const ParticleID& pid) = 0;
    virtual bool remove_voxel(const coordinate_type& coord) = 0;

    virtual bool can_move(const coordinate_type& src, const coordinate_type& dest) const = 0;
    virtual bool move(
        const coordinate_type& src, const coordinate_type& dest,
        const std::size_t candidate=0) = 0;
    virtual std::pair<coordinate_type, bool> move_to_neighbor(
        VoxelPool* const& from, VoxelPool* const& loc,
        coordinate_id_pair_type& info, const Integer nrand) = 0;

#ifdef WITH_HDF5
    virtual void save_hdf5(H5::Group* root) const = 0;
    virtual void load_hdf5(const H5::Group& root) = 0;
#endif

protected:

    virtual Integer count_voxels(const boost::shared_ptr<VoxelPool>& vp) const = 0;
    void push_voxels(std::vector<std::pair<ParticleID, Voxel> >& voxels,
            const boost::shared_ptr<MoleculePool>& voxel_pool,
            const Species& species) const;

protected:

    Real t_;
    Real voxel_radius_;

    voxel_pool_map_type voxel_pools_;
    molecule_pool_map_type molecule_pools_;

};

} // ecell4

#endif /* ECELL4_VOXELSPACEBASE_HPP */
