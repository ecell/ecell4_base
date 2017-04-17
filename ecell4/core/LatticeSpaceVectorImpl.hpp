#ifndef ECELL4_LATTICE_SPACE_VECTOR_IMPL_HPP
#define ECELL4_LATTICE_SPACE_VECTOR_IMPL_HPP

#include "LatticeSpaceBase.hpp"

namespace ecell4 {

class LatticeSpaceVectorImpl
    : public LatticeSpaceBase
{
public:

    typedef LatticeSpaceBase base_type;

    typedef base_type::coordinate_id_pair_type coordinate_id_pair_type;
    typedef base_type::coordinate_type coordinate_type;

    typedef std::vector<VoxelPool*> voxel_container;

protected:

    typedef utils::get_mapper_mf<
        Species, boost::shared_ptr<VoxelPool> >::type voxel_pool_map_type;
    typedef utils::get_mapper_mf<
        Species, boost::shared_ptr<MoleculePool> >::type molecule_pool_map_type;

public:

    LatticeSpaceVectorImpl(
        const Real3& edge_lengths, const Real& voxel_radius,
        const bool is_periodic = true);
    ~LatticeSpaceVectorImpl();

    /*
     * Space APIs
     *
     * using ParticleID, Species and Posision3
     */

    Integer num_species() const;

    virtual Real get_value(const Species& sp) const
    {
        return static_cast<Real>(num_molecules(sp));
    }

    virtual Real get_value_exact(const Species& sp) const
    {
        return static_cast<Real>(num_molecules_exact(sp));
    }

    bool has_species(const Species& sp) const;
    // bool has_species_exact(const Species& sp) const;
    virtual bool has_voxel(const ParticleID& pid) const;

    virtual bool remove_voxel(const ParticleID& pid);
    virtual bool remove_voxel(const coordinate_type& coord);

    bool update_structure(const Particle& p);

    /*
     * for Simulator
     *
     * using Species and coordinate_type
     */
    std::vector<std::pair<ParticleID, Voxel> >
        list_voxels() const;
    std::vector<std::pair<ParticleID, Voxel> >
        list_voxels(const Species& sp) const;
    std::vector<std::pair<ParticleID, Voxel> >
        list_voxels_exact(const Species& sp) const;

    virtual std::pair<ParticleID, Voxel> get_voxel(const ParticleID& pid) const;
    virtual std::pair<ParticleID, Voxel> get_voxel_at(const coordinate_type& coord) const;

    virtual Integer num_voxels_exact(const Species& sp) const;
    virtual Integer num_voxels(const Species& sp) const;
    virtual Integer num_voxels() const;

    virtual Integer num_molecules(const Species& sp) const; //XXX:

    // virtual void update_voxel(const Voxel& v);
    virtual bool update_voxel(const ParticleID& pid, const Voxel& v);

    bool add_voxels(const Species species, std::vector<std::pair<ParticleID, coordinate_type> > voxels);

    std::vector<Species> list_species() const;
    const Species& find_species(std::string name) const;
    std::vector<coordinate_type> list_coords(const Species& sp) const;
    std::vector<coordinate_type> list_coords_exact(const Species& sp) const;

    virtual bool has_molecule_pool(const Species& sp) const
    {
        molecule_pool_map_type::const_iterator itr(molecule_pools_.find(sp));
        return (itr != molecule_pools_.end());
    }

    virtual VoxelPool* find_voxel_pool(const Species& sp);
    virtual const VoxelPool* find_voxel_pool(const Species& sp) const;
    virtual MoleculePool* find_molecule_pool(const Species& sp);
    virtual const MoleculePool* find_molecule_pool(const Species& sp) const;
    // VoxelPool* find_voxel_pool(const std::string name);
    virtual VoxelPool* get_voxel_pool_at(const coordinate_type& coord) const;

    // bool update_molecule(coordinate_type coord, const Species& species);
    // bool add_molecule(const Species& sp, coordinate_type coord, const ParticleID& pid);
    virtual bool move(
        const coordinate_type& src, const coordinate_type& dest,
        const std::size_t candidate=0);
    virtual bool can_move(const coordinate_type& src, const coordinate_type& dest) const;

    std::pair<coordinate_type, bool> move_to_neighbor(
        coordinate_type coord, Integer nrand);
    std::pair<coordinate_type, bool> move_to_neighbor(
        coordinate_id_pair_type& info, Integer nrand);
    std::pair<coordinate_type, bool> move_to_neighbor(
        VoxelPool* const& from_vp, VoxelPool* const& loc,
        coordinate_id_pair_type& info, const Integer nrand);

    coordinate_type get_neighbor_boundary(
        const coordinate_type& coord, const Integer& nrand) const
    {
        coordinate_type const dest = get_neighbor(coord, nrand);
        VoxelPool* dest_vp(voxels_.at(dest));
        return (dest_vp != periodic_ ? dest : periodic_transpose(dest));
    }

    inline bool is_periodic() const
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
        // load_lattice_space(root, this, "LatticeSpaceVectorImpl");
    }
#endif

    void reset(const Real3& edge_lengths, const Real& voxel_radius,
        const bool is_periodic)
    {
        base_type::reset(edge_lengths, voxel_radius, is_periodic);

        is_periodic_ = is_periodic;
        initialize_voxels(is_periodic_);
    }

    virtual const Particle particle_at(const coordinate_type& coord) const;

    coordinate_type apply_boundary_(
        const coordinate_type& coord) const
    {
        return periodic_transpose(coord);
    }

    virtual bool make_structure_type(const Species& sp,
        Shape::dimension_kind dimension, const std::string loc);
    virtual bool make_interface_type(const Species& sp,
        Shape::dimension_kind dimension, const std::string loc);
    bool make_molecular_type(const Species& sp,
        Real radius, Real D, const std::string loc);

    virtual bool on_structure(const Voxel& v);

protected:

    VoxelPool* get_voxel_pool(const Voxel& v);

    void initialize_voxels(const bool is_periodic);

    std::pair<coordinate_type, bool> move_(
            coordinate_type from, coordinate_type to,
            const std::size_t candidate=0);
    std::pair<coordinate_type, bool> move_(
            coordinate_id_pair_type& info, coordinate_type to);
    coordinate_type get_coord(const ParticleID& pid) const;

    Integer count_voxels(const boost::shared_ptr<VoxelPool>& vp) const;

protected:

    bool is_periodic_;

    voxel_pool_map_type voxel_pools_;
    molecule_pool_map_type molecule_pools_;
    voxel_container voxels_;

    VoxelPool* vacant_;
    VoxelPool* border_;
    VoxelPool* periodic_;
};

} // ecell4

#endif
