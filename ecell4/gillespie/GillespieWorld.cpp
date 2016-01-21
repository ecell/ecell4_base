#include <vector>
#include <string>
#include <sstream>
#include <iostream>

#include <ecell4/core/SerialIDGenerator.hpp>
#include "GillespieWorld.hpp"


// using namespace std;

namespace ecell4
{

namespace gillespie
{

void GillespieWorld::set_t(const Real& t)
{
    this->cs_->set_t(t);
}

const Real GillespieWorld::t() const
{
    return this->cs_->t();
}

std::vector<Species> GillespieWorld::list_species() const
{
    return this->cs_->list_species();
}

Integer GillespieWorld::num_molecules(const Species& sp) const
{
    return this->cs_->num_molecules(sp);
}

Integer GillespieWorld::num_molecules_exact(const Species& sp) const
{
    return this->cs_->num_molecules_exact(sp);
}

Real GillespieWorld::get_value(const Species& sp) const
{
    return this->cs_->get_value(sp);
}

Real GillespieWorld::get_value_exact(const Species& sp) const
{
    return this->cs_->get_value_exact(sp);
}

void GillespieWorld::set_value(const Species& sp, const Real value)
{
    this->cs_->set_value(sp, value);
}

void GillespieWorld::add_molecules(const Species& sp, const Integer& num)
{
    this->cs_->add_molecules(sp, num);
}

void GillespieWorld::remove_molecules(const Species& sp, const Integer& num)
{
    this->cs_->remove_molecules(sp, num);
}

bool GillespieWorld::has_species(const Species& sp) const
{
    return this->cs_->has_species(sp);
}

std::vector<std::pair<ParticleID, Particle> >
    GillespieWorld::list_particles() const
{
    SerialIDGenerator<ParticleID> pidgen;
    const std::vector<Species> species_list(list_species());
    const Real3 lengths(edge_lengths());

    std::vector<std::pair<ParticleID, Particle> > retval;
    for (std::vector<Species>::const_iterator i(species_list.begin());
        i != species_list.end(); ++i)
    {
        const Integer num(num_molecules_exact(*i));

        for (Integer k(0); k < num; ++k)
        {
            const Real3 pos(
                rng_->uniform(0, lengths[0]),
                rng_->uniform(0, lengths[1]),
                rng_->uniform(0, lengths[2]));
            retval.push_back(
                std::make_pair(pidgen(), Particle(*i, pos, 0.0, 0.0)));
        }
    }
    return retval;
}

std::vector<std::pair<ParticleID, Particle> >
    GillespieWorld::list_particles_exact(const Species& sp) const
{
    SerialIDGenerator<ParticleID> pidgen;
    const Real3 lengths(edge_lengths());

    std::vector<std::pair<ParticleID, Particle> > retval;
    const Integer num(num_molecules_exact(sp));

    for (Integer k(0); k < num; ++k)
    {
        const Real3 pos(
            rng_->uniform(0, lengths[0]),
            rng_->uniform(0, lengths[1]),
            rng_->uniform(0, lengths[2]));
        retval.push_back(
            std::make_pair(pidgen(), Particle(sp, pos, 0.0, 0.0)));
    }
    return retval;
}

std::vector<std::pair<ParticleID, Particle> >
    GillespieWorld::list_particles(const Species& sp) const
{
    SerialIDGenerator<ParticleID> pidgen;
    const std::vector<Species> species_list(list_species());
    const Real3 lengths(edge_lengths());

    std::vector<std::pair<ParticleID, Particle> > retval;
    for (std::vector<Species>::const_iterator i(species_list.begin());
        i != species_list.end(); ++i)
    {
        const Integer coef(sp.count(*i));
        if (coef == 0)
        {
            continue;
        }

        const Integer num(coef * num_molecules_exact(*i));

        for (Integer k(0); k < num; ++k)
        {
            const Real3 pos(
                rng_->uniform(0, lengths[0]),
                rng_->uniform(0, lengths[1]),
                rng_->uniform(0, lengths[2]));
            retval.push_back(
                std::make_pair(pidgen(), Particle(*i, pos, 0.0, 0.0)));
        }
    }
    return retval;
}

} // gillespie

} // ecell4
