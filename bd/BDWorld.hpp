#ifndef __BD_WORLD_HPP
#define __BD_WORLD_HPP

#include <boost/scoped_ptr.hpp>

#include "../core/ParticleSpace.hpp"


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

    // // see core/Space.hpp
    // Real const& t() const
    // {
    //     return t_;
    // }

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

protected:

    boost::scoped_ptr<ParticleSpace> ps_;
};

} // bd

} // ecell4

#endif /* __BD_WORLD_HPP */
