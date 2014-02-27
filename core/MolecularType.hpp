#ifndef __ECELL4_MOLECULAR_TYPE_HPP
#define __ECELL4_MOLECULAR_TYPE_HPP

#include "MolecularTypeBase.hpp"
#include "SParticle.hpp"

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

    std::vector<SParticle> sparticles() const
    {
        std::vector<SParticle> retval;
        for (container_type::const_iterator itr(begin());
                itr != end(); ++itr)
        {
            retval.push_back(SParticle((*itr).first, &species_));
        }
        return retval;
    }

    bool is_vacant() const
    {
        return false;
    }

};

} // ecell4

#endif /* __ECELL4_MOLECULAR_TYPE_HPP */
