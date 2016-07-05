#ifndef __ECELL4_STRUCTURE_TYPE_HPP
#define __ECELL4_STRUCTURE_TYPE_HPP

#include "VoxelPool.hpp"


namespace ecell4
{

class StructureType
    : public VoxelPool
{
public:

    typedef VoxelPool base_type;
    typedef base_type::coordinate_id_pair_type coordinate_id_pair_type;
    typedef base_type::coordinate_type coordinate_type;
    typedef base_type::voxel_type_type voxel_type_type;

public:

    StructureType(
        const Species& species, VoxelPool* location,
        const Real& radius = 0.0, const Shape::dimension_kind& dimension=Shape::UNDEF)
        : base_type(species, location, radius, 0),
          dimension_(std::min(dimension, location->get_dimension()))
    {
        ;
    }

    virtual ~StructureType()
    {
        ;
    }

    virtual voxel_type_type const voxel_type() const
    {
        return STRUCTURE;
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
};

} //ecell4

#endif /* __ECELL4_STRUCTURE_TYPE_HPP */
