#include "BDContainer2D.hpp"
#include <ecell4/core/Context.hpp>

namespace ecell4
{

namespace bd
{

ParticleContainer2D::ParticleContainer2D()
{
    // do nothing
}

Integer ParticleContainer2D::num_particles() const
{
    return particles_.size();
}

Integer ParticleContainer2D::num_particles(const Species& sp) const
{
    Integer retval = 0;
    SpeciesExpressionMatcher sexp(sp);
    for(per_species_particle_id_set::const_iterator iter = particle_pool_.begin();
        iter != particle_pool_.end(); ++iter)
    {
        const Species target((*iter).first);
        if(sexp.match(target))
        {
            retval += (*iter).second.size();
        }
    }
    return retval;
}

Integer ParticleContainer2D::num_particles_exact(const Species& sp) const
{
    per_species_particle_id_set::const_iterator i =
        particle_pool_.find(sp.serial());
    if (i == particle_pool_.end())
    {
        return 0;
    }
    return (*i).second.size();
}

Integer ParticleContainer2D::num_molecules(const Species& sp) const
{
    Integer retval(0);
    SpeciesExpressionMatcher sexp(sp);
    for (per_species_particle_id_set::const_iterator i(particle_pool_.begin());
        i != particle_pool_.end(); ++i)
    {
        const Species tgt((*i).first);
        retval += sexp.count(tgt) * (*i).second.size();
    }
    return retval;
}

Integer ParticleContainer2D::num_molecules_exact(const Species& sp) const
{
    return num_particles_exact(sp);
}

std::vector<std::pair<ParticleID, Particle>>
ParticleContainer2D::list_particles() const
{
    return particles_;
}

std::vector<std::pair<ParticleID, Particle>>
ParticleContainer2D::list_particles(const Species& sp) const
{
    std::vector<std::pair<ParticleID, Particle> > retval;
    SpeciesExpressionMatcher sexp(sp);

    for (particle_container_type::const_iterator i(particles_.begin());
         i != particles_.end(); ++i)
    {
        if (sexp.match((*i).second.species()))
        {
            retval.push_back(*i);
        }
    }
    return retval;
}

std::vector<std::pair<ParticleID, Particle>>
ParticleContainer2D::list_particles_exact(const Species& sp) const
{
    std::vector<std::pair<ParticleID, Particle> > retval;

    for (particle_container_type::const_iterator i(particles_.begin());
         i != particles_.end(); ++i)
    {
        if ((*i).second.species() == sp)
        {
            retval.push_back(*i);
        }
    }
    return retval;
}



}// bd
}// ecell4
