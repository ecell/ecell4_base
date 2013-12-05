#ifndef __ECELL4_MOLECULAR_TYPE_BASE_HPP
#define __ECELL4_MOLECULAR_TYPE_BASE_HPP

#include <vector>
#include "Species.hpp"
#include "Identifier.hpp"

namespace ecell4
{

class MolecularTypeBase
{
public:
    typedef std::vector<std::pair<Integer, ParticleID> > container_type;

public:
    virtual void addVoxel(Integer coord, const ParticleID pid)
    {
        throw "addVoxel(Integer, ParticleID) is not supported.";
    }

    virtual bool removeVoxel(Integer coord)
    {
        throw "removeVoxel(Integer) is not supported.";
    }

    virtual const Species& species() const
    {
        throw "species() is not supported.";
    }

    virtual const container_type& voxels() const
    {
        throw "voxels() is not supported.";
    }

};

} // ecell4

#endif
