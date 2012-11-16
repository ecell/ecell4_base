#include <cmath>
// #include <gsl/gsl_pow_int.h>

#include "exceptions.hpp"
#include "ParticleSpace.hpp"


namespace ecell4
{

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

bool ParticleSpaceVectorImpl::update_particle(
    ParticleID const& pid, Particle const& p)
{
    index_map_type::iterator i(index_map_.find(pid));
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
    index_map_type::iterator i(index_map_.find(pid));
    if (i == index_map_.end())
    {
        // throw not_found("Particle not found");
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
        throw not_found("Particle not found");
    }

    return particles_[(*i).second];
}

Real ParticleSpaceVectorImpl::distance_sq(
    Position3 const& p1, Position3 const& p2) const
{
    // not implemented yet
    return 0.0;
}

std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
ParticleSpaceVectorImpl::get_particles_within_radius(
    Position3 const& pos, Real const& radius) const
{
    // Real const rsq(gsl_pow_2(radius));
    Real const rsq(radius * radius);
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

} // ecell4
