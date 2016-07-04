#ifndef __ECELL4_VACANT_TYPE_HPP
#define __ECELL4_VACANT_TYPE_HPP

#include "MolecularTypeBase.hpp"

namespace ecell4
{

class VacantType
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

    ~VacantType()
    {
        ; // do nothing
    }

    virtual voxel_type_type const voxel_type() const
    {
        return VACANT;
    }

    bool with_voxels() const
    {
        return false;
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

#endif /* __ECELL4_VACANT_TYPE_HPP */
