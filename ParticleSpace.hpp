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

    virtual Position3 const& edge_lengths() const = 0;

    virtual Integer num_species() const = 0;
    virtual Integer num_particles() const = 0;
    virtual Integer num_particles(Species const& species) const = 0;

    virtual bool update_particle(ParticleID const& pid, Particle const& p) = 0;
    virtual bool remove_particle(ParticleID const& pid) = 0;

    virtual std::pair<ParticleID, Particle>
    get_particle(ParticleID const& pid) const = 0;
    virtual std::vector<std::pair<ParticleID, Particle> >
    get_particles(Species const& species) const = 0;

    virtual Position3 apply_boundary(Position3 const& pos) const = 0;
    virtual Real distance_sq(
        Position3 const& pos1, Position3 const& pos2) const = 0;

    inline Real distance(Position3 const& pos1, Position3 const& pos2) const
    {
        return std::sqrt(distance_sq(pos1, pos2));
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
    typedef typename container_type::size_type index_type;
    typedef std::map<ParticleID, index_type> index_map_type;

    ParticleSpaceVectorImpl(Position3 const& edge_lengths)
        : edge_lengths_(edge_lengths)
    {
        ;
    }

    Position3 const& edge_lengths() const
    {
        return edge_lengths_;
    }

    Integer num_species() const;
    Integer num_particles() const;
    Integer num_particles(Species const& species) const;

    bool update_particle(ParticleID const& pid, Particle const& p);
    bool remove_particle(ParticleID const& pid);

    std::pair<ParticleID, Particle> get_particle(ParticleID const& pid) const;
    std::vector<std::pair<ParticleID, Particle> >
    get_particles(Species const& species) const;

    Position3 apply_boundary(Position3 const& pos) const = 0;
    Real distance_sq(Position3 const& pos1, Position3 const& pos2) const;

    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
    get_particles_within_radius(
        Position3 const& pos, Real const& radius) const;

protected:

    Position3 edge_lengths_;
    container_type particles_;
    index_map_type index_map_;
};

} // ecell4

#endif /* __PARTICLE_SPACE_HPP */
