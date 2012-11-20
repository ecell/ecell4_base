#include <cmath>

#include <ecell4/core/exceptions.hpp>
#include "BDSimulator.hpp"


namespace ecell4
{

namespace bd
{

Real I_bd_3d(Real const& sigma, Real const& t, Real const& D)
{
    const Real sqrtPi(std::sqrt(M_PI));

    const Real Dt(D * t);
    const Real Dt2(Dt + Dt);
    const Real sqrtDt(std::sqrt(Dt));
    const Real sigmasq(sigma * sigma);

    const Real term1(1 / (3 * sqrtPi));
    const Real term2(sigmasq - Dt2);
    const Real term3(Dt2 - 3 * sigmasq);
    const Real term4(sqrtPi * sigmasq * sigma * erfc(sigma / sqrtDt));

    const Real result(
        term1 * (-sqrtDt * (term2 * std::exp(-sigmasq / Dt) + term3) + term4));
    return result;
}

bool BDSimulator::attempt_reaction(
    ParticleID const& pid, Particle const& particle)
{
    ReactionRuleVector reaction_rules(
        (*model_).query_reaction_rules(particle.species()));
    if (reaction_rules.size() == 0)
    {
        return false;
    }

    Real const rnd((*state_).rng.uniform(0, 1));
    Real prob(0);
    for (ReactionRuleVector::const_iterator i(reaction_rules.begin());
         i != reaction_rules.end(); ++i)
    {
        ReactionRule const& rr(*i);
        prob += rr.k() * dt();
        if (prob > rnd)
        {
            SpeciesVector const& products(rr.products());
            switch (products.size())
            {
            case 0:
                break;
            case 1:
                break;
            case 2:
                break;
            default:
                throw NotImplemented(
                    "more than two products are not allowed");
                break;
            }
            return true;
        }
    }

    return false;
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

    Real const D1(particle1.D()), D2(particle2.D());
    Real const r01(particle1.radius() + particle2.radius());
    Real const rnd((*state_).rng.uniform(0, 1));
    Real prob(0);

    for (ReactionRuleVector::const_iterator i(reaction_rules.begin());
         i != reaction_rules.end(); ++i)
    {
        ReactionRule const& rr(*i);
        prob += rr.k() * dt() / (
            (I_bd_3d(r01, dt(), D1) + I_bd_3d(r01, dt(), D2)) * 4 * M_PI);

        if (prob >= 1)
        {
            throw std::runtime_error(
                "the total reaction probability exceeds 1."
                " the step interval is too long");
        }
        if (prob > rnd)
        {
            SpeciesVector const& products(rr.products());
            switch (products.size())
            {
            case 0:
                break;
            case 1:
                break;
            case 2:
                break;
            default:
                throw NotImplemented(
                    "more than two products are not allowed");
                break;
            }
            return true;
        }
    }

    return false;
}

Position3 BDSimulator::draw_displacement_3d(Particle const& particle)
{
    Real const sigma(std::sqrt(2 * particle.D() * dt()));
    return Position3(
        (*state_).rng.gaussian(0, sigma),
        (*state_).rng.gaussian(0, sigma),
        (*state_).rng.gaussian(0, sigma));
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
