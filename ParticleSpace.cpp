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

typename ParticleSpace::particle_id_pair_type
ParticleSpaceVectorImpl::new_particle(Species const& sp, Position const& pos)
{
    return std::make_pair(ParticleID(), Particle());
}

bool ParticleSpaceVectorImpl::update_particle(
    particle_id_pair_type const& pidpair)
{
    return true;
}

bool ParticleSpaceVectorImpl::remove_particle(ParticleID const& pid)
{
    return true;
}

} // ecell4
