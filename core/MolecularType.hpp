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
    typedef MolecularTypeBase::container_type container_type;

public:
    MolecularType(const std::string& name = "")
        : species_(name)
    {
    }

    MolecularType(const Species& species) : species_(species)
    {
    }

    void addVoxel(Integer coord, const ParticleID pid);
    bool removeVoxel(Integer coord);
    const Species& species() const;
    std::vector<SParticle> sparticles() const;
    const container_type& voxels() const;

    container_type::iterator begin()
    {
        return voxels_.begin();
    }

    container_type::iterator end()
    {
        return voxels_.end();
    }

    container_type::const_iterator begin() const
    {
        return voxels_.begin();
    }

    container_type::const_iterator end() const
    {
        return voxels_.end();
    }

    container_type::iterator find(Integer coord)
    {
        container_type::iterator itr;
        for (itr = voxels_.begin(); itr != voxels_.end(); ++itr)
        {
            if ((*itr).first == coord)
            {
                break;
            }
        }
        return itr;
    }

protected:
    Species species_;
    container_type voxels_;

};

} // ecell4

#endif /* __ECELL4_MOLECULAR_TYPE_HPP */
