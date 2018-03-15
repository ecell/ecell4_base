#ifndef ECELL4_VACANT_TYPE_HPP
#define ECELL4_VACANT_TYPE_HPP
#include <boost/shared_ptr.hpp>
#include "VoxelPool.hpp"

namespace ecell4
{

class VacantType
    : public VoxelPool
{
public:

    ~VacantType()
    {
        ; // do nothing
    }

    virtual voxel_type_type const voxel_type() const
    {
        return VACANT;
    }

    static
    boost::shared_ptr<VacantType>
    allocate(const Shape::dimension_kind& dimension=Shape::THREE)
    {
        return boost::shared_ptr<VacantType>(new VacantType(dimension));
    }

    const Shape::dimension_kind get_dimension() const
    {
        return dimension_;
    }

private:

    const Shape::dimension_kind dimension_;

private:

    typedef VoxelPool base_type;

    VacantType(const Shape::dimension_kind& dimension)
        : base_type(Species("VACANT", "0", "0"), boost::weak_ptr<VoxelPool>(), 0, 0),
          dimension_(dimension)
    {
        ; // do nothing
    }
};

} // ecell4

#endif /* ECELL4_VACANT_TYPE_HPP */
