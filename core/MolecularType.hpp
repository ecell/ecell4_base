#ifndef __ECELL4_MOLECULAR_TYPE_HPP
#define __ECELL4_MOLECULAR_TYPE_HPP

#include "MolecularTypeBase.hpp"

namespace ecell4
{

class MolecularType
    : public MolecularTypeBase
{

public:
    typedef MolecularTypeBase::particle_info particle_info;
    typedef MolecularTypeBase::container_type container_type;

public:
    MolecularType(const std::string& name = "")
        : MolecularTypeBase(Species(name))
    {
    }

    MolecularType(const Species& species)
        : MolecularTypeBase(species)
    {
    }

    ~MolecularType()
    {
    }

    bool is_vacant() const
    {
        return false;
    }

};

} // ecell4

#endif /* __ECELL4_MOLECULAR_TYPE_HPP */
