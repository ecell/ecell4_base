#ifndef __ECELL4_MOLECULAR_TYPE_BASE_HPP
#define __ECELL4_MOLECULAR_TYPE_BASE_HPP

#include <vector>
#include "Species.hpp"
#include "Shape.hpp"
#include "Identifier.hpp"
#include "RandomNumberGenerator.hpp"
#include "Identifier.hpp"
// #include "LatticeSpace.hpp"


namespace ecell4
{

class MolecularTypeBase
{
public:

    typedef Integer coordinate_type;
    typedef std::pair<coordinate_type, ParticleID> coord_id_pair;

    typedef std::vector<coord_id_pair> container_type;
    typedef container_type::const_iterator const_iterator;
    typedef container_type::iterator iterator;

public:

    MolecularTypeBase(
        const Species& species, MolecularTypeBase* location,
        const Real& radius, const Real& D)
        : species_(species), location_(location),
        radius_(radius), D_(D), dimension_(Shape::UNDEF)
    {
        ;
    }

    virtual ~MolecularTypeBase()
    {
        ;
    }

    virtual bool is_vacant() const = 0;

    virtual bool with_voxels() const
    {
        return true;
    }

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

    void set_dimension(const Shape::dimension_kind dimension)
    {
        dimension_ = dimension;
    }

    virtual const Shape::dimension_kind get_dimension() const
    {
        return dimension_;
    }

    virtual void add_voxel_without_checking(const coord_id_pair& info)
    {
        voxels_.push_back(info);
    }

    virtual void replace_voxel(
        const coordinate_type& from_coord,
        const coord_id_pair& to_info)
    {
        container_type::iterator itr(find(from_coord));
        if (itr == voxels_.end())
        {
            throw NotFound("no corresponding coordinate was found.");
        }

        (*itr) = to_info;
    }

    virtual void replace_voxel(
        const coordinate_type& from_coord,
        const coordinate_type& to_coord)
    {
        container_type::iterator itr(find(from_coord));
        if (itr == voxels_.end())
        {
            throw NotFound("no corresponding coordinate was found.");
        }

        (*itr).first = to_coord;
    }

    virtual coord_id_pair pop(const coordinate_type& coord)
    {
        container_type::iterator position(this->find(coord));
        const coord_id_pair info(*position);
        this->remove_voxel(position);
        return info;
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

    void remove_voxel(const container_type::iterator& position)
    {
        // voxels_.erase(position);
        (*position) = voxels_.back();
        voxels_.pop_back();
    }

    void swap(
        const container_type::iterator& a, const container_type::iterator& b)
    {
        if (a == b)
        {
            return;
        }

        const container_type::value_type info(*b);
        (*b) = (*a);
        (*a) = info;
    }

    coord_id_pair& at(const Integer& index)
    {
        return voxels_.at(index);
    }

    coord_id_pair const& at(const Integer& index) const
    {
        return voxels_.at(index);
    }

    coord_id_pair& operator[](const Integer& n)
    {
        return voxels_[n];
    }

    coord_id_pair const& operator[](const Integer& n) const
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

    const ParticleID find_particle_id(const coordinate_type& coord) const
    {
        container_type::const_iterator i(this->find(coord));
        if (i == voxels_.end())
        {
            throw NotFound("No corresponding ParticleID was found.");
        }
        return (*i).second;
    }

    container_type::iterator find(const ParticleID& pid)
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

    container_type::const_iterator find(const ParticleID& pid) const
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

    container_type::iterator find(coordinate_type coord)
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

    container_type::const_iterator find(coordinate_type coord) const
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

protected:

    const Species species_;
    MolecularTypeBase* location_;
    Real radius_, D_;
    Shape::dimension_kind dimension_;

    container_type voxels_;
};

} // ecell4

#endif
