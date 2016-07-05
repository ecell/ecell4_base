#ifndef __ECELL4_INTERFACE_TYPE_HPP
#define __ECELL4_INTERFACE_TYPE_HPP

#include "MolecularTypeBase.hpp"


namespace ecell4
{

class InterfaceType
    : public VoxelPool
{
public:

    typedef VoxelPool base_type;
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

private:

    const Shape::dimension_kind dimension_;

    // const Species species_;
    // VoxelPool* location_;
    // Real radius_, D_;
    // container_type voxels_;  // typedef std::vector<coordinate_id_pair_type> container_type;
};

} // ecell4

#endif /* __ECELL4_INTERFACE_TYPE_HPP */
