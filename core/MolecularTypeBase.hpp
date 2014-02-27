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
    typedef std::pair<Coord, ParticleID> particle_info;
    typedef std::vector<particle_info> container_type;

public:
    MolecularTypeBase(const Species& species) : species_(species)
    {
    }

    virtual ~MolecularTypeBase()
    {
    }

    virtual bool is_vacant() const = 0;

    const Species& species() const
    {
        return species_;
    }

    void addVoxel(particle_info info)
    {
        container_type::iterator itr(find(info.first));
        if (itr != voxels_.end())
        {
            voxels_.erase(itr);
        }
        voxels_.push_back(info);
    }

    bool removeVoxel(Coord coord)
    {
        container_type::iterator itr(find(coord));
        if (itr != voxels_.end())
        {
            voxels_.erase(itr);
            return true;
        }
        return false;
    }

    const container_type& voxels() const
    {
        return voxels_;
    }

    container_type& voxels()
    {
        return voxels_;
    }

    container_type::iterator begin()
    {
        return voxels_.begin();
    }

    container_type::const_iterator begin() const
    {
        return voxels_.begin();
    }

    container_type::iterator end()
    {
        return voxels_.end();
    }

    container_type::const_iterator end() const
    {
        return voxels_.end();
    }

    container_type::iterator find(Coord coord)
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

    container_type::const_iterator find(Coord coord) const
    {
        container_type::const_iterator itr;
        for (itr = voxels_.begin(); itr != voxels_.end(); ++itr)
        {
            if ((*itr).first == coord)
            {
                break;
            }
        }
        return itr;
    }

    container_type::iterator find(ParticleID pid)
    {
        container_type::iterator itr;
        for (itr = voxels_.begin(); itr != voxels_.end(); ++itr)
        {
            if ((*itr).second == pid)
            {
                break;
            }
        }
        return itr;
    }

    container_type::const_iterator find(ParticleID pid) const
    {
        container_type::const_iterator itr;
        for (itr = voxels_.begin(); itr != voxels_.end(); ++itr)
        {
            if ((*itr).second == pid)
            {
                break;
            }
        }
        return itr;
    }

protected:
    const Species species_;
    container_type voxels_;

};

} // ecell4

#endif
