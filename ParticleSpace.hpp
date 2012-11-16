#ifndef __PARTICLE_SPACE_HPP
#define __PARTICLE_SPACE_HPP

#include "types.hpp"
#include "Position3.hpp"
#include "Particle.hpp"
#include "Species.hpp"
#include "Space.hpp"


namespace ecell4
{

class ParticleSpace
    : public Space
{
public:

    virtual Integer num_species() const = 0;

    virtual Integer num_particles() const = 0;
    virtual bool update_particle(ParticleID const& pid, Particle const& p) = 0;
    virtual bool remove_particle(ParticleID const& pid) = 0;

    virtual std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
    get_particles_within_radius(
        Position3 const& pos, Real const& radius) const = 0;
};

class ParticleSpaceVectorImpl
    : public ParticleSpace
{
public:

    ParticleSpaceVectorImpl()
    {
        ;
    }

    Integer num_species() const;

    Integer num_particles() const;
    bool update_particle(ParticleID const& pid, Particle const& p);
    bool remove_particle(ParticleID const& pid);

    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
    get_particles_within_radius(
        Position3 const& pos, Real const& radius) const;
};

} // ecell4

#endif /* __PARTICLE_SPACE_HPP */
