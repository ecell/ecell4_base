#ifndef __ECELL4_ODE_ODE_WORLD_HPP
#define __ECELL4_ODE_ODE_WORLD_HPP

#include <ecell4/core/Species.hpp>

namespace ecell4
{

namespace ode
{

class ODEWorld
{
protected:

    typedef std::vector<Real> num_molecules_container_type;
    typedef std::vector<Species> species_container_type;
    typedef utils::get_mapper_mf<
        Species, num_molecules_container_type::size_type>::type species_map_type;

public:

    ODEWorld(const Real& volume)
        : volume_(volume), t_(0.0)
    {
        ;
    }

    // SpaceTraits

    const Real& t() const
    {
        return t_;
    }

    void set_t(const Real& t)
    {
        if (t < 0.0)
        {
            throw std::invalid_argument("the time must be positive.");
        }
        t_ = t;
    }

    // CompartmentSpaceTraits

    const Real& volume() const
    {
        return volume_;
    }

    Integer num_species(void)
    {
        return static_cast<Integer>(species_.size());
    }

    bool has_species(const Species& sp)
    {
        species_map_type::const_iterator i(index_map_.find(sp));
        return (i != index_map_.end());
    }

    Real num_molecules(const Species& sp)
    {
        species_map_type::const_iterator i(index_map_.find(sp));
        if (i == index_map_.end())
        {
            throw NotFound("Species not found");
        }

        return num_molecules_[(*i).second];
    }

    // CompartmentSpace member functions

    void set_volume(const Real& volume)
    {
        if (volume <= 0.0)
        {
            throw std::invalid_argument("The volume must be positive.");
        }

        volume_ = volume;
    }

    void add_species(const Species& sp)
    {
        species_map_type::const_iterator i(index_map_.find(sp));
        if (i != index_map_.end())
        {
            throw AlreadyExists("Species already exists");
        }

        index_map_.insert(std::make_pair(sp, num_molecules_.size()));
        species_.push_back(sp);
        num_molecules_.push_back(0);
    }

    void remove_species(const Species& sp)
    {
        species_map_type::iterator i(index_map_.find(sp));
        if (i == index_map_.end())
        {
            throw NotFound("Species not found");
        }

        species_map_type::mapped_type
            idx((*i).second), last_idx(num_molecules_.size() - 1);
        if (idx != last_idx)
        {
            const species_container_type::size_type
                idx_(static_cast<species_container_type::size_type>(idx)),
                last_idx_(
                    static_cast<species_container_type::size_type>(last_idx));
            const Species& last_sp(species_[last_idx_]);
            species_[idx_] = last_sp;
            num_molecules_[idx] = num_molecules_[last_idx];
            index_map_[last_sp] = idx;
        }

        species_.pop_back();
        num_molecules_.pop_back();
        index_map_.erase(sp);
    }

    void add_molecules(const Species& sp, const Real& num)
    {
        set_num_molecules(sp, num_molecules(sp) + num);
    }

    void remove_molecules(const Species& sp, const Real& num)
    {
        set_num_molecules(sp, num_molecules(sp) - num);
    }

    // Optional members

    void set_num_molecules(const Species& sp, const Real& num)
    {
        if (num < 0)
        {
            throw std::invalid_argument(
                "The number of molecules must be positive.");
        }

        species_map_type::const_iterator i(index_map_.find(sp));
        if (i == index_map_.end())
        {
            throw NotFound("Species not found");
        }

        num_molecules_[(*i).second] = num;
    }

protected:

    Real volume_;
    Real t_;

    num_molecules_container_type num_molecules_;
    species_container_type species_;
    species_map_type index_map_;
};

} // ode

} // ecell4

#endif /* __ECELL4_ODE_ODE_WORLD_HPP */
