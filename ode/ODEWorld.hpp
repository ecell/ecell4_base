#ifndef __ODEWORLD_HPP
#define __ODEWORLD_HPP

#include <ecell4/core/Species.hpp>

namespace ecell4
{

namespace ode
{

class ODEWorld
{
public:

    ODEWorld(Real const& volume)
        : volume_(volume)
    {
        ;
    }

    Real const& volume() const
    {
        return volume_;
    }

    void set_volume(Real const& volume)
    {
        if (volume <= 0)
        {
            throw std::invalid_argument("The volume must be positive.");
        }

        volume_ = volume;
    }

    Integer num_species(void)
    {
        return static_cast<Integer>(species_.size());
    }

    bool has_species(Species const &sp)
    {
        species_map_type::const_iterator i(index_map_.find(sp));
        return (i != index_map_.end());
    }

    void add_species(Species const &sp)
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

    void remove_species(Species const &sp)
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
            species_container_type::size_type const
                idx_(static_cast<species_container_type::size_type>(idx)),
                last_idx_(static_cast<species_container_type::size_type>(last_idx));
            Species const& last_sp(species_[last_idx_]);
            species_[idx_] = last_sp;
            num_molecules_[idx] = num_molecules_[last_idx];
            index_map_[last_sp] = idx;
        }

        species_.pop_back();
        num_molecules_.pop_back();
        index_map_.erase(sp);
    }

    Real num_molecules(Species const& sp)
    {
        species_map_type::const_iterator i(index_map_.find(sp));
        if (i == index_map_.end())
        {
            throw NotFound("Species not found");
        }

        return num_molecules_[(*i).second];
    }

    void set_num_molecules(Species const& sp, Real const& num)
    {
        if (num < 0)
        {
            throw std::invalid_argument("The number of molecules must be positive.");
        }

        species_map_type::const_iterator i(index_map_.find(sp));
        if (i == index_map_.end())
        {
            throw NotFound("Species not found");
        }

        num_molecules_[(*i).second] = num;
    }

protected:

    typedef std::vector<Real> num_molecules_container_type;
    typedef std::vector<Species> species_container_type;

    typedef utils::get_mapper_mf<
        Species, num_molecules_container_type::size_type>::type species_map_type;

    Real volume_;

    num_molecules_container_type num_molecules_;
    species_container_type species_;
    species_map_type index_map_;
};

} // ode

} // ecell4

#endif // __ODEWORLD_HPP
