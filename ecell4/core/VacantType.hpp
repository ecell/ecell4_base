#ifndef ECELL4_VACANT_TYPE_HPP
#define ECELL4_VACANT_TYPE_HPP

#include "VoxelPool.hpp"

namespace ecell4
{

class VacantType
    : public VoxelPool
{
public:

    typedef VoxelPool base_type;
    typedef base_type::coordinate_id_pair_type coordinate_id_pair_type;
    typedef base_type::coordinate_type coordinate_type;
    typedef base_type::voxel_type_type voxel_type_type;

public:

    ~VacantType()
    {
        ; // do nothing
    }

    virtual voxel_type_type const voxel_type() const
    {
        return VACANT;
    }

    static VacantType& getInstance()
    {
        static VacantType instance;
        return instance;
    }

    const Shape::dimension_kind get_dimension() const
    {
        return Shape::THREE;
    }

private:

    VacantType()
        : base_type(Species("VACANT", "0", "0"), NULL, 0, 0)
    {
        ; // do nothing
    }
};

} // ecell4

#endif /* ECELL4_VACANT_TYPE_HPP */
