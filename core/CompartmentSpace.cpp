#include <stdexcept>

#include "exceptions.hpp"
#include "CompartmentSpace.hpp"


namespace ecell4
{

Real const& CompartmentSpaceVectorImpl::volume() const
{
    return volume_;
}

void CompartmentSpaceVectorImpl::set_volume(Real volume)
{
    if (volume <= 0)
    {
        throw std::invalid_argument("The volume must be positive.");
    }

    volume_ = volume;
}

void CompartmentSpaceVectorImpl::add_species(Species const& sp)
{
    index_map_type::const_iterator i(index_map_.find(sp));
    if (i != index_map_.end())
    {
        throw AlreadyExists("Species already exists");
    }

    index_map_.insert(std::make_pair(sp, num_molecules_.size()));
    species_.push_back(sp);
    num_molecules_.push_back(0);
}

void CompartmentSpaceVectorImpl::remove_species(Species const& sp)
{
    index_map_type::iterator i(index_map_.find(sp));
    if (i == index_map_.end())
    {
        throw NotFound("Species not found");
    }

    index_type idx((*i).second), last_idx(num_molecules_.size() - 1);
    if (idx != last_idx)
    {
        Species const& last_sp(species_[last_idx]);
        species_[idx] = last_sp;
        num_molecules_[idx] = num_molecules_[last_idx];
        index_map_[last_sp] = idx;
    }

    species_.pop_back();
    num_molecules_.pop_back();
    index_map_.erase(sp);
}

bool CompartmentSpaceVectorImpl::has_species(Species const& sp) const
{
    index_map_type::const_iterator i(index_map_.find(sp));
    return (i != index_map_.end());
}

Integer CompartmentSpaceVectorImpl::num_species() const
{
    return static_cast<Integer>(species_.size());
}

Integer CompartmentSpaceVectorImpl::num_molecules(Species const& sp) const
{
    index_map_type::const_iterator i(index_map_.find(sp));
    if (i == index_map_.end())
    {
        throw NotFound("Species not found");
    }

    return num_molecules_[(*i).second];
}

void CompartmentSpaceVectorImpl::add_molecules(
    Species const& sp, Integer const& num)
{
    if (num <= 0)
    {
        throw std::invalid_argument("The number of molecules must be positive.");
    }

    index_map_type::const_iterator i(index_map_.find(sp));
    if (i == index_map_.end())
    {
        throw NotFound("Species not found");
    }

    num_molecules_[(*i).second] += num;
}

void CompartmentSpaceVectorImpl::remove_molecules(
    Species const& sp, Integer const& num)
{
    if (num <= 0)
    {
        throw std::invalid_argument("The number of molecules must be positive.");
    }

    index_map_type::const_iterator i(index_map_.find(sp));
    if (i == index_map_.end())
    {
        throw NotFound("Species not found");
    }

    if (num_molecules_[(*i).second] < num)
    {
        throw std::invalid_argument(
            "The number of molecules cannot be negative.");
    }

    num_molecules_[(*i).second] -= num;
}

} // ecell4
