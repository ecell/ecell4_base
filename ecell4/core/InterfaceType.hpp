#ifndef ECELL4_INTERFACE_TYPE_HPP
#define ECELL4_INTERFACE_TYPE_HPP

#include "VoxelPool.hpp"


namespace ecell4
{

class InterfaceType
    : public VoxelPool
{
public:

    typedef VoxelPool base_type;
    typedef base_type::coordinate_id_pair_type coordinate_id_pair_type;
    typedef base_type::coordinate_type coordinate_type;
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

    const Shape::dimension_kind get_dimension() const
    {
        return dimension_;
    }

private:

    const Shape::dimension_kind dimension_;
};

} // ecell4

#endif /* ECELL4_INTERFACE_TYPE_HPP */
