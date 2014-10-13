#ifndef __ECELL4_MOLECULAR_TYPE_HPP
#define __ECELL4_MOLECULAR_TYPE_HPP

#include "MolecularTypeBase.hpp"
#include "VacantType.hpp"

namespace ecell4
{

class MolecularType
    : public MolecularTypeBase
{
public:

    typedef MolecularTypeBase base_type;
    typedef base_type::particle_info particle_info;
    typedef base_type::private_coordinate_type private_coordinate_type;
    typedef base_type::container_type container_type;

public:

    MolecularType(const std::string& name = "")
        : MolecularTypeBase(Species(name), &(VacantType::getInstance()), 0, 0)
    {
        ;
    }

    MolecularType(const Species& species, const Real& radius = 0.0,
            const Real& D = 0.0)
        : MolecularTypeBase(species, &(VacantType::getInstance()), radius, D)
    {
        ;
    }

    MolecularType(const Species& species, MolecularTypeBase* location,
            const Real& radius = 0.0, const Real& D = 0.0)
        : MolecularTypeBase(species, location, radius, D)
    {
        ;
    }

    ~MolecularType()
    {
        ;
    }

    bool is_vacant() const
    {
        return false;
    }
};

} // ecell4

#endif /* __ECELL4_MOLECULAR_TYPE_HPP */
