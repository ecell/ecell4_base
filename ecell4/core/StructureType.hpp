#ifndef __ECELL4_STRUCTURE_TYPE_HPP
#define __ECELL4_STRUCTURE_TYPE_HPP

#include "MolecularTypeBase.hpp"
#include "VacantType.hpp"

namespace ecell4
{

class StructureType
    : public MolecularTypeBase
{

public:

    StructureType(const Species& species, const Real& radius = 0.0)
        : MolecularTypeBase(species, &(VacantType::getInstance()), radius, 0)
    {
    }

    ~StructureType()
    {
    }

    void addVoxel(LatticeSpace::private_coordinate_type coord)
    {
        MolecularTypeBase::addVoxel(particle_info(coord, ParticleID()));
    }

    bool is_vacant() const
    {
        return false;
    }

};

} //ecell4

#endif /* __ECELL4_STRUCTURE_TYPE_HPP */
