#ifndef __ECELL4_MOLECULAR_TYPE_BASE_HPP
#define __ECELL4_MOLECULAR_TYPE_BASE_HPP

#include <vector>
#include "Species.hpp"
#include "Identifier.hpp"
#include "RandomNumberGenerator.hpp"
#include "LatticeSpace.hpp"

namespace ecell4
{

class MolecularTypeBase
{
public:

    typedef LatticeSpace::particle_info particle_info;
    typedef std::vector<particle_info> container_type;

public:

    MolecularTypeBase(const Species& species, MolecularTypeBase* location,
            const Real& radius, const Real& D)
        : species_(species), location_(location), radius_(radius), D_(D)
    {
        ;
    }

    virtual ~MolecularTypeBase()
    {
        ;
    }

    virtual bool is_vacant() const = 0;

    const Species& species() const
    {
        return species_;
    }

    MolecularTypeBase* location() const
    {
        return location_;
    }

    Real& radius()
    {
        return radius_;
    }

    const Real& radius() const
    {
        return radius_;
    }

    Real& D()
    {
        return D_;
    }

    const Real& D() const
    {
        return D_;
    }

    virtual void addVoxel(particle_info info)
    {
        container_type::iterator itr(find(info.first));
        if (itr != voxels_.end())
        {
            voxels_.erase(itr);
        }
        voxels_.push_back(info);
    }

    virtual bool removeVoxel(LatticeSpace::private_coordinate_type coord)
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

    particle_info& at(const Integer& index)
    {
        return voxels_.at(index);
    }

    const Integer size() const
    {
        return voxels_.size();
    }

    void shuffle(RandomNumberGenerator& rng)
    {
        ecell4::shuffle(rng, voxels_);
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

    container_type::iterator find(LatticeSpace::private_coordinate_type coord)
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

    container_type::const_iterator find(LatticeSpace::private_coordinate_type coord) const
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
    MolecularTypeBase* location_;
    Real radius_, D_;

    container_type voxels_;
};

} // ecell4

#endif
