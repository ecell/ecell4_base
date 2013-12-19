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
    typedef std::pair<Coord, ParticleID> particle_info;
    typedef std::vector<particle_info> container_type;

public:
    MolecularType(const std::string& name = "")
        : species_(name)
    {
    }

    MolecularType(const Species& species) : species_(species)
    {
    }

    ~MolecularType()
    {
    }

    void addVoxel(particle_info info);
    bool removeVoxel(Coord coord);
    const Species& species() const;
    std::vector<SParticle> sparticles() const;
    const container_type& voxels() const;

    bool is_vacant() const
    {
        return false;
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
    Species species_;
    container_type voxels_;

};

} // ecell4

#endif /* __ECELL4_MOLECULAR_TYPE_HPP */
