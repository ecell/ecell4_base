#ifndef __ECELL4_INTERFACE_TYPE_HPP
#define __ECELL4_INTERFACE_TYPE_HPP

#include "MolecularTypeBase.hpp"


namespace ecell4
{

class InterfaceType
    : public MolecularType
    // : public MolecularTypeBase
{
public:

    typedef MolecularType base_type;
    // typedef MolecularTypeBase base_type;
    typedef base_type::coordinate_id_pair_type coordinate_id_pair_type;
    typedef base_type::coordinate_type coordinate_type;
    typedef base_type::container_type container_type;
    typedef base_type::iterator iterator;
    typedef base_type::const_iterator const_iterator;
    typedef base_type::voxel_type_type voxel_type_type;

public:

    InterfaceType(
        const Species& sp, VoxelPool* location,
        const Real& radius = 0, const Shape::dimension_kind& dimension=Shape::UNDEF)
        : base_type(sp, location, radius, 0),
        dimension_(std::min(dimension, location->get_dimension()))
    {
        ;
    }

    ~InterfaceType()
    {
        ; // do nothing
    }

    voxel_type_type const voxel_type() const
    {
        return INTERFACE;
    }

    bool with_voxels() const
    {
        return false;
    }

    const Shape::dimension_kind get_dimension() const
    {
        return dimension_;
    }

    virtual void add_voxel_without_checking(const coordinate_id_pair_type& info)
    {
        ; // do nothing
    }

    virtual void replace_voxel(
        const coordinate_type& from_coord,
        const coordinate_id_pair_type& to_info)
    {
        ; // do nothing
    }

    virtual void replace_voxel(
        const coordinate_type& from_coord,
        const coordinate_type& to_coord,
        const std::size_t candidate=0)
    {
        ; // do nothing
    }

    virtual coordinate_id_pair_type pop(const coordinate_type& coord)
    {
        // container_type::iterator position(this->find(coord));
        // const coordinate_id_pair_type info(*position);
        // this->remove_voxel(position);
        // return info;
        return coordinate_id_pair_type(ParticleID(), coord);
    }

    bool remove_voxel_if_exists(const coordinate_type& coord)
    {
        return true; // just return true
    }

    virtual const ParticleID get_particle_id(const coordinate_type& coord) const
    {
        return ParticleID();
    }

private:

    const Shape::dimension_kind dimension_;

    // const Species species_;
    // VoxelPool* location_;
    // Real radius_, D_;
    // container_type voxels_;  // typedef std::vector<coordinate_id_pair_type> container_type;
};

} // ecell4

#endif /* __ECELL4_INTERFACE_TYPE_HPP */
