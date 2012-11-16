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
    species_.erase(std::unique(species_.begin(), species_.end()),
                   species_.end());
    return static_cast<Integer>(species_.size());
}

Integer ParticleSpaceVectorImpl::num_particles() const
{
    return static_cast<Integer>(particles_.size());
}

bool ParticleSpaceVectorImpl::update_particle(
    ParticleID const& pid, Particle const& p)
{
    return true;
}

bool ParticleSpaceVectorImpl::remove_particle(ParticleID const& pid)
{
    return true;
}

Real ParticleSpaceVectorImpl::distance_sq(
    Position3 const& p1, Position3 const& p2) const
{
    return 0.0;
}

std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
ParticleSpaceVectorImpl::get_particles_within_radius(
    Position3 const& pos, Real const& radius) const
{
    return std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >();
}

} // ecell4
