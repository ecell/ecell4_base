#ifndef __ECELL4_VACANT_TYPE_HPP
#define __ECELL4_VACANT_TYPE_HPP

#include "MolecularTypeBase.hpp"

namespace ecell4
{

class VacantType
    : public MolecularTypeBase
{

public:
    ~VacantType()
    {
    }

    static VacantType& getInstance()
    {
        static VacantType instance;
        return instance;
    }

    void addVoxel(particle_info info)
    {
    }

    bool removeVoxel(LatticeSpace::private_coordinate_type coord)
    {
        return true;
    }

    bool is_vacant() const
    {
        return true;
    }
private:
    VacantType() : MolecularTypeBase(Species("VACANT", "0"), NULL, 0, 0)
    {
    }

};

} // ecell4

#endif /* __ECELL4_VACANT_TYPE_HPP */
