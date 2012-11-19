#include <cmath>
#include <stdexcept>
// #include <gsl/gsl_pow_int.h>

#include "Exceptions.hpp"
#include "ParticleSpace.hpp"


namespace ecell4
{

Real pow_2(Real const& a)
{
    // return gsl_pow_2(a);
    return a * a;
}

Integer ParticleSpaceVectorImpl::num_species() const
{
    std::vector<Species> species_;
    for (container_type::const_iterator i(particles_.begin());
         i != particles_.end(); ++i)
    {
        species_.push_back((*i).second.species());
    }

    std::sort(species_.begin(), species_.end());
    return static_cast<Integer>(
        species_.end() - std::unique(species_.begin(), species_.end()));
}

Integer ParticleSpaceVectorImpl::num_particles() const
{
    return static_cast<Integer>(particles_.size());
}

Integer ParticleSpaceVectorImpl::num_particles(Species const& species) const
{
    return static_cast<Integer>(get_particles(species).size());
}

bool ParticleSpaceVectorImpl::has_particle(ParticleID const& pid) const
{
    index_map_type::const_iterator i(index_map_.find(pid));
    return (i != index_map_.end());
}

bool ParticleSpaceVectorImpl::update_particle(
    ParticleID const& pid, Particle const& p)
{
    index_map_type::const_iterator i(index_map_.find(pid));
    if (i == index_map_.end())
    {
        particles_[(*i).second] = std::make_pair(pid, p);
        return true;
    }
    else
    {
        index_map_[pid] = particles_.size();
        particles_.push_back(std::make_pair(pid, p));
        return false;
    }
}

bool ParticleSpaceVectorImpl::remove_particle(ParticleID const& pid)
{
    index_map_type::const_iterator i(index_map_.find(pid));
    if (i == index_map_.end())
    {
        // throw NotFound("Particle not found");
        return false;
    }

    index_type idx((*i).second), last_idx(particles_.size() - 1);
    if (idx != last_idx)
    {
        std::pair<ParticleID, Particle> const& last(particles_[last_idx]);
        particles_[idx] = last;
        index_map_[last.first] = idx;
    }

    particles_.pop_back();
    index_map_.erase((*i).first);
    return true;
}

std::pair<ParticleID, Particle> ParticleSpaceVectorImpl::get_particle(
    ParticleID const& pid) const
{
    index_map_type::const_iterator i(index_map_.find(pid));
    if (i == index_map_.end())
    {
        throw NotFound("Particle not found");
    }

    return particles_[(*i).second];
}

std::vector<std::pair<ParticleID, Particle> >
ParticleSpaceVectorImpl::get_particles() const
{
    return particles_;
}

std::vector<std::pair<ParticleID, Particle> >
ParticleSpaceVectorImpl::get_particles(Species const& species) const
{
    std::vector<std::pair<ParticleID, Particle> > retval;

    for (container_type::const_iterator i(particles_.begin());
         i != particles_.end(); ++i)
    {
        if ((*i).second.species() == species)
        {
            retval.push_back(*i);
        }
    }

    return retval;
}

Position3 ParticleSpaceVectorImpl::apply_boundary(Position3 const& pos) const
{
    return modulo(pos, edge_lengths_);
}

Real ParticleSpaceVectorImpl::distance_sq(
    Position3 const& p1, Position3 const& p2) const
{
    Real retval(0);
    for (Position3::size_type dim(0); dim < 3; ++dim)
    {
        const Real edge_length(edge_lengths_[dim]);
        const Real diff(p1[dim] - p2[dim]), half(edge_length * 0.5);

        if (diff > half)
        {
            retval += pow_2(diff - edge_length);
        }
        else if (diff < -half)
        {
            retval += pow_2(diff + edge_length);
        }
        else
        {
            retval += pow_2(diff);
        }
    }
    return retval;
}

std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
ParticleSpaceVectorImpl::get_particles_within_radius(
    Position3 const& pos, Real const& radius) const
{
    Real const rsq(pow_2(radius));
    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> > retval;

    for (container_type::const_iterator i(particles_.begin());
         i != particles_.end(); ++i)
    {
        Real const dsq(distance_sq((*i).second.position(), pos));
        if (dsq <= rsq)
        {
            retval.push_back(std::make_pair(*i, std::sqrt(dsq)));
        }
    }

    return retval;
}

void ParticleSpaceVectorImpl::set_edge_lengths(Position3 const& edge_lengths)
{
    for (Position3::size_type dim(0); dim < 3; ++dim)
    {
        if (edge_lengths[dim] <= 0)
        {
            throw std::invalid_argument("the edge length must be positive.");
        }
    }

    edge_lengths_ = edge_lengths;
}

} // ecell4
