#ifndef __ECELL4_OFFLATTICE_SPACE_HPP
#define __ECELL4_OFFLATTICE_SPACE_HPP

#include "LatticeSpace.hpp"

namespace ecell4 {

class OffLatticeSpace
    : public LatticeSpace
{
public:

    typedef LatticeSpace super;
    typedef super::coordinate_type coordinate_type;
    typedef super::coordinate_id_pair_type coordinate_id_pair_type;

public:

    OffLatticeSpace(const Real& voxel_radius)
        : t_(0.0), voxel_radius_(voxel_radius)
    {
        ;
    }

    virtual ~OffLatticeSpace()
    {
        ;
    }

    virtual Real3 actual_lengths() const = 0;

    virtual std::vector<Species> list_species() const = 0;

    virtual Integer num_voxels_exact(const Species& sp) const = 0;
    virtual Integer num_voxels(const Species& sp) const = 0;
    virtual Integer num_voxels() const = 0;
    virtual bool has_voxel(const ParticleID& pid) const = 0;

    virtual std::vector<std::pair<ParticleID, Voxel> > list_voxels() const = 0;
    virtual std::vector<std::pair<ParticleID, Voxel> > list_voxels(const Species& sp) const = 0;
    virtual std::vector<std::pair<ParticleID, Voxel> > list_voxels_exact(const Species& sp) const = 0;

    virtual bool update_voxel(const ParticleID& pid, const Voxel& v) = 0;
    virtual bool update_voxel_without_checking(const ParticleID& pid, const Voxel& v)
    {
        throw NotSupported(
            "update_voxel_without_chekcing(const ParticleID&, const Voxel&) is not supported by this space class");
    }

    virtual std::pair<ParticleID, Voxel> get_voxel(const ParticleID& pid) const = 0;
    virtual std::pair<ParticleID, Voxel> get_voxel(const coordinate_type& coord) const = 0;
    virtual bool remove_voxel(const ParticleID& pid) = 0;
    virtual bool remove_voxel(const coordinate_type& coord) = 0;
    virtual bool move(
        const coordinate_type& src, const coordinate_type& dest,
        const std::size_t candidate=0) = 0;

    // Not pure virtual function
    virtual bool can_move(const coordinate_type& src, const coordinate_type& dest) const;

    virtual const Particle particle_at(const coordinate_type& coord) const = 0;

    virtual bool has_molecule_pool(const Species& sp) const = 0;
    virtual VoxelPool* find_voxel_pool(const Species& sp) = 0;
    virtual const VoxelPool* find_voxel_pool(const Species& sp) const = 0;
    virtual MoleculePool* find_molecule_pool(const Species& sp) = 0;
    virtual const MoleculePool* find_molecule_pool(const Species& sp) const = 0;
    virtual VoxelPool* find_voxel_pool(const coordinate_type& coord) const = 0;

    // Not pure virtual function
    virtual bool make_structure_type(const Species& sp,
        Shape::dimension_kind dimension, const std::string loc);
    virtual bool make_interface_type(const Species& sp,
        Shape::dimension_kind dimension, const std::string loc);

    virtual bool on_structure(const Voxel& v) = 0;

    virtual std::pair<coordinate_type, bool> move_to_neighbor(
        VoxelPool* const& from_mt, VoxelPool* const& loc,
        coordinate_id_pair_type& info, const Integer nrand) = 0;

    virtual const Integer col_size() const = 0;
    virtual const Integer row_size() const = 0;
    virtual const Integer layer_size() const = 0;

}

}

#endif /* __ECELL4_OFFLATTICE_SPACE_HPP */
