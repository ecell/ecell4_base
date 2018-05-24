#ifndef ECELL4_STRUCTURE_TYPE_HPP
#define ECELL4_STRUCTURE_TYPE_HPP

#include "VoxelPool.hpp"


namespace ecell4
{

static inline Shape::dimension_kind
calc_dimension(boost::weak_ptr<VoxelPool> location, const Shape::dimension_kind& dimension)
{
    if (location.expired())
        return dimension;

    return std::min(dimension, location.lock()->get_dimension());
}

class StructureType
    : public VoxelPool
{
private:

    typedef VoxelPool base_type;

public:

    StructureType(
        const Species& species, boost::weak_ptr<VoxelPool> location,
        const Real& radius = 0.0, const Shape::dimension_kind& dimension=Shape::UNDEF)
        : base_type(species, location, radius, 0),
          dimension_(calc_dimension(location, dimension)),
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

    const Shape::dimension_kind get_dimension() const
    {
        return dimension_;
    }

private:

    const Shape::dimension_kind dimension_;
};

} //ecell4

#endif /* ECELL4_STRUCTURE_TYPE_HPP */
