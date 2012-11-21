#ifndef __BD_WORLD_HPP
#define __BD_WORLD_HPP

#include <boost/scoped_ptr.hpp>

#include <ecell4/core/SerialIDGenerator.hpp>
#include <ecell4/core/ParticleSpace.hpp>


namespace ecell4
{

namespace bd
{

class BDWorld
{
public:

    BDWorld(Position3 const& edge_lengths)
        : ps_(new ParticleSpaceVectorImpl(edge_lengths))
    {
        ;
    }

    ParticleID new_particle_id()
    {
        return pidgen_();
    }

    Real const& t() const
    {
        return (*ps_).t();
    }

    void set_t(Real const& t)
    {
        (*ps_).set_t(t);
    }

    Position3 const& edge_lengths() const
    {
        return (*ps_).edge_lengths();
    }

    Integer num_species() const
    {
        return (*ps_).num_species();
    }

    Integer num_particles() const
    {
        return (*ps_).num_particles();
    }

    Integer num_particles(Species const& species) const
    {
        return (*ps_).num_particles(species);
    }

    bool has_particle(ParticleID const& pid) const
    {
        return (*ps_).has_particle(pid);
    }

    bool update_particle(ParticleID const& pid, Particle const& p)
    {
        return (*ps_).update_particle(pid, p);
    }

    bool remove_particle(ParticleID const& pid)
    {
        return (*ps_).remove_particle(pid);
    }

    std::pair<ParticleID, Particle>
    get_particle(ParticleID const& pid) const
    {
        return (*ps_).get_particle(pid);
    }

    std::vector<std::pair<ParticleID, Particle> > get_particles() const
    {
        return (*ps_).get_particles();
    }

    std::vector<std::pair<ParticleID, Particle> >
    get_particles(Species const& species) const
    {
        return (*ps_).get_particles(species);
    }

    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
    get_particles_within_radius(
        Position3 const& pos, Real const& radius) const
    {
        return (*ps_).get_particles_within_radius(pos, radius);
    }

    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
    get_particles_within_radius(
        Position3 const& pos, Real const& radius, ParticleID const& ignore) const
    {
        return (*ps_).get_particles_within_radius(pos, radius, ignore);
    }

    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
    get_particles_within_radius(
        Position3 const& pos, Real const& radius,
        ParticleID const& ignore1, ParticleID const& ignore2) const
    {
        return (*ps_).get_particles_within_radius(pos, radius, ignore1, ignore2);
    }

    inline Position3 periodic_transpose(
        Position3 const& pos1, Position3 const& pos2) const
    {
        return (*ps_).periodic_transpose(pos1, pos2);
    }

    inline Position3 apply_boundary(Position3 const& pos) const
    {
        return (*ps_).apply_boundary(pos);
    }

    inline Real distance_sq(Position3 const& pos1, Position3 const& pos2) const
    {
        return (*ps_).distance_sq(pos1, pos2);
    }

    inline Real distance(Position3 const& pos1, Position3 const& pos2) const
    {
        return (*ps_).distance(pos1, pos2);
    }

protected:

    boost::scoped_ptr<ParticleSpace> ps_;
    SerialIDGenerator<ParticleID> pidgen_;
};

} // bd

} // ecell4

#endif /* __BD_WORLD_HPP */
