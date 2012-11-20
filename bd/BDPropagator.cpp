#include <ecell4/core/exceptions.hpp>

#include "BDPropagator.hpp"


namespace ecell4
{

namespace bd
{

Position3 random_displacement_3d(
    RandomNumberGenerator& rng, Real const& t, Real const& D)
{
    Real const sigma(std::sqrt(2 * D * t));
    return Position3(
        rng.gaussian(0, sigma), rng.gaussian(0, sigma), rng.gaussian(0, sigma));
}

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

bool BDPropagator::operator()()
{
    if (queue_.empty())
    {
        return false;
    }

    ParticleID const pid(queue_.back().first);
    queue_.pop_back();
    Particle particle(world_.get_particle(pid).second);

    if (attempt_reaction(pid, particle))
    {
        return true;
    }

    Real const D(particle.D());
    if (D == 0)
    {
        return true;
    }

    Position3 const newpos(
        particle.position() + draw_displacement(particle));
    Particle particle_to_update(
        particle.species(), newpos, particle.radius(), particle.D());
    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
        overlapped(world_.get_particles_within_radius(
                       newpos, particle.radius(), pid));

    switch (overlapped.size())
    {
    case 0:
        world_.update_particle(pid, particle_to_update);
        return true;
    case 1:
    {
        std::pair<ParticleID, Particle> closest(
            (*(overlapped.begin())).first);
        if (attempt_reaction(
                pid, particle_to_update, closest.first, closest.second))
        {
            return true;
        }
        return true;
    }
    default:
        return true;
    }
}

bool BDPropagator::attempt_reaction(
    ParticleID const& pid, Particle const& particle)
{
    ReactionRuleVector reaction_rules(
        model_.query_reaction_rules(particle.species()));
    if (reaction_rules.size() == 0)
    {
        return false;
    }

    Real const rnd(rng().uniform(0, 1));
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

bool BDPropagator::attempt_reaction(
    ParticleID const& pid1, Particle const& particle1,
    ParticleID const& pid2, Particle const& particle2)
{
    ReactionRuleVector reaction_rules(
        model_.query_reaction_rules(
            particle1.species(), particle2.species()));
    if (reaction_rules.size() == 0)
    {
        return false;
    }

    Real const D1(particle1.D()), D2(particle2.D());
    Real const r01(particle1.radius() + particle2.radius());
    Real const rnd(rng().uniform(0, 1));
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

} // bd

} // ecell4
