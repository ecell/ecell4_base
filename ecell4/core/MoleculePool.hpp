#ifndef ECELL4_MOLECULE_POOL_HPP
#define ECELL4_MOLECULE_POOL_HPP

#include "VoxelPool.hpp"

namespace ecell4 {

class MoleculePool
    : public VoxelPool
{
public:

    typedef VoxelPool base_type;

    typedef std::vector<coordinate_id_pair_type> container_type;
    typedef container_type::const_iterator const_iterator;
    typedef container_type::iterator iterator;

public:

    MoleculePool(
        const Species& species, boost::weak_ptr<VoxelPool> location,
        const Real& radius=0.0, const Real& D=0.0)
        : base_type(species, location, radius, D)
    {
        ;
    }

    virtual ~MoleculePool()
    {
        ;
    }

    virtual voxel_type_type const voxel_type() const
    {
        return DEFAULT;
    }

public:

    virtual void add_voxel(const coordinate_id_pair_type& info)
    {
        voxels_.push_back(info);
    }

    virtual void replace_voxel(
        const coordinate_type& from_coord, const coordinate_type& to_coord,
        const std::size_t candidate = 0)
    {
        container_type::iterator itr(find(from_coord, candidate));
        if (itr == voxels_.end())
        {
            std::cerr << "from_coord = " << from_coord << std::endl;
            throw NotFound("no corresponding coordinate was found.");
        }

        (*itr).coordinate = to_coord;
    }

    virtual bool remove_voxel_if_exists(const coordinate_type& coord)
    {
        container_type::iterator itr(find(coord));
        if (itr != voxels_.end())
        {
            this->remove_voxel(itr);
            return true;
        }
        return false;
    }

    virtual const ParticleID get_particle_id(const coordinate_type& coord) const
    {
        container_type::const_iterator i(this->find(coord));
        if (i == voxels_.end())
        {
            throw NotFound("No corresponding ParticleID was found.");
        }
        return (*i).pid;
    }

public:

    void remove_voxel(const container_type::iterator& position)
    {
        // voxels_.erase(position);
        (*position) = voxels_.back();
        voxels_.pop_back();
    }

    coordinate_id_pair_type pop(const coordinate_type& coord)
    {
        container_type::iterator position(this->find(coord));
        const coordinate_id_pair_type info(*position);
        this->remove_voxel(position);
        return info;
    }

    void replace_voxel(
        const coordinate_type& from_coord, const coordinate_id_pair_type& to_info)
    {
        container_type::iterator itr(find(from_coord));
        if (itr == voxels_.end())
        {
            throw NotFound("no corresponding coordinate was found.");
        }

        (*itr) = to_info;
    }

    void swap(const container_type::iterator& a, const container_type::iterator& b)
    {
        if (a == b)
        {
            return;
        }

        const container_type::value_type info(*b);
        (*b) = (*a);
        (*a) = info;
    }

    coordinate_id_pair_type& at(const Integer& index)
    {
        return voxels_.at(index);
    }

    coordinate_id_pair_type const& at(const Integer& index) const
    {
        return voxels_.at(index);
    }

    coordinate_id_pair_type& operator[](const Integer& n)
    {
        return voxels_[n];
    }

    coordinate_id_pair_type const& operator[](const Integer& n) const
    {
        return voxels_[n];
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

    container_type::iterator find(const ParticleID& pid)
    {
        container_type::iterator itr;
        for (itr = voxels_.begin(); itr != voxels_.end(); ++itr)
        {
            if ((*itr).pid == pid)
            {
                break;
            }
        }
        return itr;
    }

    container_type::const_iterator find(const ParticleID& pid) const
    {
        container_type::const_iterator itr;
        for (itr = voxels_.begin(); itr != voxels_.end(); ++itr)
        {
            if ((*itr).pid == pid)
            {
                break;
            }
        }
        return itr;
    }

protected:

    container_type::iterator find(
        coordinate_type coord, const std::size_t candidate = 0)
    {
        container_type::iterator itr;
        if (candidate < voxels_.size())
        {
            itr = voxels_.begin() + candidate;
            if ((*itr).coordinate == coord)
                return itr;
        }
        for (itr = voxels_.begin(); itr != voxels_.end(); ++itr)
        {
            if ((*itr).coordinate == coord)
            {
                break;
            }
        }
        return itr;
    }

    container_type::const_iterator find(
        coordinate_type coord, const std::size_t candidate = 0) const
    {
        container_type::const_iterator itr;
        if (candidate < voxels_.size())
        {
            itr = voxels_.begin() + candidate;
            if ((*itr).coordinate == coord)
                return itr;
        }
        for (itr = voxels_.begin(); itr != voxels_.end(); ++itr)
        {
            if ((*itr).coordinate == coord)
            {
                break;
            }
        }
        return itr;
    }

protected:

    container_type voxels_;
};

} // ecell4

#endif
