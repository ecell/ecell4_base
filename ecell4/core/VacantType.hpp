#ifndef ECELL4_VACANT_TYPE_HPP
#define ECELL4_VACANT_TYPE_HPP
#include <boost/shared_ptr.hpp>
#include "StructureType.hpp"

namespace ecell4
{

class VacantType
    : public StructureType
{
private:

    typedef StructureType base_type;

    VacantType()
        : base_type(Species("", "0", "0"), boost::weak_ptr<VoxelPool>(), 0)
    {
        ; // do nothing
    }

public:

    ~VacantType()
    {
        ; // do nothing
    }

    voxel_type_type const voxel_type() const
    {
        return VACANT;
    }

    static
    boost::shared_ptr<VacantType>
    allocate()
    {
        return boost::shared_ptr<VacantType>(new VacantType());
    }
};

} // ecell4

#endif /* ECELL4_VACANT_TYPE_HPP */
