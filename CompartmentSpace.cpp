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
    if (volume > 0)
    {
        volume_ = volume;
    }
}

void CompartmentSpaceVectorImpl::add_species(Species const& sp)
{
    index_map_type::const_iterator i(index_map_.find(sp));
    if (i != index_map_.end())
    {
        index_map_.insert(std::make_pair(sp, num_of_molecules_.size()));
        species_.push_back(sp);
        num_of_molecules_.push_back(0);
    }
}

Integer CompartmentSpaceVectorImpl::num_of_molecules(Species const& sp)
{
    index_map_type::const_iterator i(index_map_.find(sp));
    if (i == index_map_.end())
    {
        throw not_found("Species not found");
    }

    return num_of_molecules_[(*i).second];
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
        throw not_found("Species not found");
    }

    num_of_molecules_[(*i).second] += num;
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
        throw not_found("Species not found");
    }

    if (num_of_molecules_[(*i).second] < num)
    {
        throw std::invalid_argument(
            "The number of molecules cannot be negative.");
    }

    num_of_molecules_[(*i).second] -= num;
}

} // ecell4
