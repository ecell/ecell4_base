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

    typedef std::pair<ParticleID, Particle> particle_id_pair_type;

    virtual Integer num_species() const = 0;

    virtual Integer num_particles() const = 0;
    virtual particle_id_pair_type new_particle(
        Species const& sp, Position3 const& pos) = 0;
    virtual bool update_particle(particle_id_pair_type const& pidpair) = 0;
    virtual bool remove_particle(ParticleID const& pid) = 0;
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
    particle_id_pair_type new_particle(Species const& sp, Position3 const& pos);
    bool update_particle(particle_id_pair_type const& pidpair);
    bool remove_particle(ParticleID const& pid);
};

} // ecell4

#endif /* __PARTICLE_SPACE_HPP */
