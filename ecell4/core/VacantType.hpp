#ifndef ECELL4_VACANT_TYPE_HPP
#define ECELL4_VACANT_TYPE_HPP
#include "StructureType.hpp"
#include <memory>

namespace ecell4
{

class VacantType : public StructureType
{
private:
    typedef StructureType base_type;

    VacantType() : base_type(Species("", 0, 0), std::weak_ptr<VoxelPool>())
    {
        ; // do nothing
    }

public:
    ~VacantType()
    {
        ; // do nothing
    }

    voxel_type_type const voxel_type() const { return VACANT; }

    static std::shared_ptr<VacantType> allocate()
    {
        return std::shared_ptr<VacantType>(new VacantType());
    }
};

} // namespace ecell4

#endif /* ECELL4_VACANT_TYPE_HPP */
