#ifndef ECELL4_INTERFACE_TYPE_HPP
#define ECELL4_INTERFACE_TYPE_HPP

#include "StructureType.hpp"


namespace ecell4
{

class InterfaceType
    : public StructureType
{
private:

    typedef StructureType base_type;

public:

    InterfaceType(
        const Species& sp, boost::weak_ptr<VoxelPool> location,
        const Real& radius = 0, const Shape::dimension_kind& dimension=Shape::UNDEF)
        : base_type(sp, location, radius, dimension)
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

};

} // ecell4

#endif /* ECELL4_INTERFACE_TYPE_HPP */
