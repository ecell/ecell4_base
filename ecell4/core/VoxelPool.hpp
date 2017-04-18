#ifndef ECELL4_MOLECULAR_TYPE_BASE_HPP
#define ECELL4_MOLECULAR_TYPE_BASE_HPP

#include <vector>
#include "Species.hpp"
#include "Shape.hpp"
#include "Identifier.hpp"
#include "RandomNumberGenerator.hpp"
#include "Voxel.hpp"


namespace ecell4
{

class VoxelPool
{
public:

    typedef Voxel::coordinate_type coordinate_type;

    // typedef std::pair<coordinate_type, ParticleID> coordinate_id_pair_type;
    typedef struct coordinate_id_pair_type
        {
            coordinate_id_pair_type(ParticleID const& pid, coordinate_type const& coordinate)
                : coordinate(coordinate), pid(pid)
            {}

            coordinate_type coordinate;
            ParticleID pid;

            bool operator==(const coordinate_id_pair_type& rhs) const
            {
                return pid == rhs.pid &&
                    coordinate == rhs.coordinate;
            }

            bool operator!=(const coordinate_id_pair_type& rhs) const
            {
                return pid != rhs.pid
                    || coordinate != rhs.coordinate;
            }

            bool operator<(const coordinate_id_pair_type& rhs) const
            {
                return coordinate < rhs.coordinate
                    || (coordinate == rhs.coordinate &&
                        pid < rhs.pid);
            }

            bool operator>=(const coordinate_id_pair_type& rhs) const
            {
                return coordinate > rhs.coordinate
                    || (coordinate == rhs.coordinate &&
                        pid >= rhs.pid);
            }

            bool operator>(const coordinate_id_pair_type& rhs) const
            {
                return coordinate > rhs.coordinate
                    || (coordinate == rhs.coordinate &&
                        pid > rhs.pid);
            }

            bool operator<=(const coordinate_id_pair_type& rhs) const
            {
                return coordinate < rhs.coordinate
                    || (coordinate == rhs.coordinate &&
                        pid <= rhs.pid);
            }
        } coordinate_id_pair_type;

public:

    typedef enum
    {
        DEFAULT,
        VACANT,
        STRUCTURE,
        INTERFACE
    } voxel_type_type;

public:

    VoxelPool(
        const Species& species, VoxelPool* location,
        const Real& radius, const Real& D)
        : species_(species), location_(location),
        radius_(radius), D_(D)
    {
        ;
    }


    virtual ~VoxelPool()
    {
        ;
    }

    virtual voxel_type_type const voxel_type() const = 0;

    virtual const Shape::dimension_kind get_dimension() const
    {
        return Shape::UNDEF;
    }

public:

    bool is_vacant() const
    {
        return voxel_type() == VACANT;
    }

    bool is_structure() const
    {
        return voxel_type() == STRUCTURE;
    }

    bool is_interface() const
    {
        return voxel_type() == INTERFACE;
    }

    const Species& species() const
    {
        return species_;
    }

    VoxelPool* location() const
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

public:

    virtual void add_voxel(const coordinate_id_pair_type& info)
    {
        if (info.pid != ParticleID())
        {
            throw NotSupported("No ParticleID is allowed.");
        }

        ; // do nothing
    }

    virtual void replace_voxel(
        const coordinate_type& from_coord, const coordinate_type& to_coord,
        const std::size_t candidate = 0)
    {
        ; // do nothing
    }

    virtual coordinate_id_pair_type pop(const coordinate_type& coord)
    {
        return coordinate_id_pair_type(ParticleID(), coord);
    }

    virtual bool remove_voxel_if_exists(const coordinate_type& coord)
    {
        return true;
    }

    virtual const ParticleID get_particle_id(const coordinate_type& coord) const
    {
        return ParticleID();
    }

protected:

    const Species species_;
    VoxelPool* location_;
    Real radius_, D_;
};

class MoleculePool
    : public VoxelPool
{
public:

    typedef VoxelPool base_type;
    typedef base_type::coordinate_type coordinate_type;
    typedef base_type::coordinate_id_pair_type coordinate_id_pair_type;
    typedef base_type::voxel_type_type voxel_type_type;

    typedef std::vector<coordinate_id_pair_type> container_type;
    typedef container_type::const_iterator const_iterator;
    typedef container_type::iterator iterator;

public:

    MoleculePool(
        const Species& species, VoxelPool* location,
        const Real& radius, const Real& D)
        : base_type(species, location, radius, D)
    {
        ;
    }

    virtual ~MoleculePool()
    {
        ;
    }

    // virtual voxel_type_type const voxel_type() const = 0;

    // virtual const Shape::dimension_kind get_dimension() const = 0;

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
