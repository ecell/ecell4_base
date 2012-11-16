#ifndef __PARTICLE_SPACE_HPP
#define __PARTICLE_SPACE_HPP

#include <cmath>

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

    virtual Real const& edge_length() const = 0;

    virtual Integer num_species() const = 0;

    virtual Integer num_particles() const = 0;
    virtual bool update_particle(ParticleID const& pid, Particle const& p) = 0;
    virtual bool remove_particle(ParticleID const& pid) = 0;

    virtual Real distance_sq(Position3 const& p1, Position3 const& p2) const = 0;

    inline Real distance(Position3 const& p1, Position3 const& p2) const
    {
        return std::sqrt(distance_sq(p1, p2));
    }

    virtual std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
    get_particles_within_radius(
        Position3 const& pos, Real const& radius) const = 0;
};

class ParticleSpaceVectorImpl
    : public ParticleSpace
{
public:

    typedef std::vector<std::pair<ParticleID, Particle> > container_type;

    ParticleSpaceVectorImpl(Real const& edge_length)
        : edge_length_(edge_length)
    {
        ;
    }

    Real const& edge_length() const
    {
        return edge_length_;
    }

    Integer num_species() const;

    Integer num_particles() const;
    bool update_particle(ParticleID const& pid, Particle const& p);
    bool remove_particle(ParticleID const& pid);

    Real distance_sq(Position3 const& p1, Position3 const& p2) const;

    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
    get_particles_within_radius(
        Position3 const& pos, Real const& radius) const;

protected:

    Real edge_length_;
    container_type particles_;
    std::map<ParticleID, typename container_type::size_type> index_map_;
};

} // ecell4

#endif /* __PARTICLE_SPACE_HPP */
