#include "VacantType.hpp"

namespace ecell4
{

void VacantType::addVoxel(Voxel *voxel)
{
    voxel->ptr_mt = this;
}

} // ecell4
