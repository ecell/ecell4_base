#ifndef __ECELL4_MOLECULAR_TYPE_BASE_HPP
#define __ECELL4_MOLECULAR_TYPE_BASE_HPP

#include <vector>
#include "Species.hpp"
#include "Identifier.hpp"

namespace ecell4
{

typedef Integer Coord;

class MolecularTypeBase
{
public:
    virtual ~MolecularTypeBase()
    {
    }

    virtual void addVoxel(std::pair<Coord, ParticleID> info)
    {
        throw "addVoxel(Coord, ParticleID) is not supported.";
    }

    virtual bool removeVoxel(Coord coord)
    {
        throw "removeVoxel(Coord) is not supported.";
    }

    virtual const Species& species() const
    {
        throw "species() const is not supported.";
    }

    virtual const std::vector<std::pair<Coord, ParticleID> >& voxels() const
    {
        throw "voxels() const is not supported.";
    }

    virtual bool is_vacant() const = 0;

    virtual std::vector<std::pair<Coord, ParticleID> >::iterator
        begin()
    {
        throw "begin() is not supported.";
    }

    virtual std::vector<std::pair<Coord, ParticleID> >::const_iterator
        begin() const
    {
        throw "begin() const is not supported.";
    }

    virtual std::vector<std::pair<Coord, ParticleID> >::iterator
        end()
    {
        throw "end() is not supported.";
    }

    virtual std::vector<std::pair<Coord, ParticleID> >::const_iterator
        end() const
    {
        throw "end() const is not supported.";
    }

    virtual std::vector<std::pair<Coord, ParticleID> >::iterator
        find(Coord coord)
    {
        throw "find(Coord) is not supported.";
    }

    virtual std::vector<std::pair<Coord, ParticleID> >::const_iterator
        find(Coord coord) const
    {
        throw "find(Coord) const is not supported.";
    }

    virtual std::vector<std::pair<Coord, ParticleID> >::iterator
        find(ParticleID pid)
    {
        throw "find(ParticleID) is not supported.";
    }

    virtual std::vector<std::pair<Coord, ParticleID> >::const_iterator
        find(ParticleID pid) const
    {
        throw "find(ParticleID) const is not supported.";
    }

};

} // ecell4

#endif
