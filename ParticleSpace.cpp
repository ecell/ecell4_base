#include "ParticleSpace.hpp"


namespace ecell4
{

Integer ParticleSpaceVectorImpl::num_species() const
{
    return 0;
}

Integer ParticleSpaceVectorImpl::num_particles() const
{
    return 0;
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

std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
ParticleSpaceVectorImpl::get_particles_within_radius(
    Position3 const& pos, Real const& radius) const
{
    return std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >();
}

} // ecell4
