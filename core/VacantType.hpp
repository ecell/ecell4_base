#ifndef __ECELL4_VACANT_TYPE_HPP
#define __ECELL4_VACANT_TYPE_HPP

#include "MolecularTypeBase.hpp"

namespace ecell4
{

class VacantType
    : public MolecularTypeBase
{

public:
    typedef MolecularTypeBase::container_type container_type;

public:
    void addVoxel(Voxel *voxel);

};

} // ecell4

#endif /* __ECELL4_VACANT_TYPE_HPP */
