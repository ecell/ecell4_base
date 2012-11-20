#include "BDSimulator.hpp"


namespace ecell4
{

namespace bd
{

bool BDSimulator::attempt_reaction(
    ParticleID const& pid, Particle const& particle)
{
    ReactionRuleVector reaction_rules(
        (*model_).query_reaction_rules(particle.species()));
    if (reaction_rules.size() == 0)
    {
        return false;
    }

    return true;
}

bool BDSimulator::attempt_reaction(
    ParticleID const& pid1, Particle const& particle1,
    ParticleID const& pid2, Particle const& particle2)
{
    ReactionRuleVector reaction_rules(
        (*model_).query_reaction_rules(
            particle1.species(), particle2.species()));
    if (reaction_rules.size() == 0)
    {
        return false;
    }

    return true;
}

Position3 BDSimulator::draw_displacement(Particle const& particle)
{
    Real const sigma(std::sqrt(2 * particle.D() * dt()));
    return Position3(
        rng_.gaussian(0, sigma),
        rng_.gaussian(0, sigma),
        rng_.gaussian(0, sigma));
}

void BDSimulator::step()
{
    Real const dt_(dt());
    std::vector<std::pair<ParticleID, Particle> >
        queue_((*world_).get_particles());
    shuffle((*state_).rng, queue_);

    while (!queue_.empty())
    {
        ParticleID const pid(queue_.back().first);
        queue_.pop_back();
        Particle particle((*world_).get_particle(pid).second);

        if (attempt_reaction(pid, particle))
        {
            continue;
        }

        Real const D(particle.D());
        if (D == 0)
        {
            continue;
        }

        Position3 const newpos(
            particle.position() + draw_displacement(particle));
        Particle particle_to_update(
            particle.species(), newpos, particle.radius(), particle.D());
        std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
            overlapped((*world_).get_particles_within_radius(
                           newpos, particle.radius(), pid));

        switch (overlapped.size())
        {
        case 0:
            (*world_).update_particle(pid, particle_to_update);
            break;
        case 1:
        {
            std::pair<ParticleID, Particle> closest(
                (*(overlapped.begin())).first);
            if (attempt_reaction(
                    pid, particle_to_update, closest.first, closest.second))
            {
                continue;
            }
            continue;
        }
        default:
            continue;
        }
    }

    set_t(t() + dt_);
    ++(*state_).num_steps;
}

bool BDSimulator::step(Real const& upto)
{
    Real const t_(t()), dt_(dt());
    Real const next_time(t_ + dt_);

    if (upto > next_time)
    {
        step();
        return false;
    }
    else
    {
        set_dt(next_time - t_);
        step();
        set_dt(dt_);
        return true;
    }
}

} // bd

} // ecell4
