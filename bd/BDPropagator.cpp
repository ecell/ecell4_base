#include <ecell4/core/exceptions.hpp>

#include "BDPropagator.hpp"


namespace ecell4
{

namespace bd
{

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
                remove_particle(pid);
                break;
            case 1:
            {
                Species const& species_new(*(products.begin()));
                // Real const radius_new(species_new.radius());
                // Real const D_new(species_new.D());
                Real const radius_new(particle.radius());
                Real const D_new(particle.D());

                std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
                    overlapped(world_.get_particles_within_radius(
                                   particle.position(), radius_new, pid));
                if (overlapped.size() > 0)
                {
                    // throw NoSpace("");
                    return false;
                }

                Particle particle_to_update(
                    species_new, particle.position(), radius_new, D_new);
                world_.update_particle(pid, particle_to_update);
                break;
            }
            case 2:
            {
                Species const& species_new1(products[0]);
                Species const& species_new2(products[1]);
                // Real const D1(species_new1.D()), D2(species_new2.D());
                // Real const radius1(species_new1.radius()),
                //     radius2(species_new2.radius());
                Real const D1(particle.D()), D2(particle.D());
                Real const radius1(particle.radius()), radius2(particle.radius());

                Real const D12(D1 + D2);
                Real const r12(radius1 + radius2);
                Position3 newpos1, newpos2;
                Integer i(max_retry_count_);
                while (true)
                {
                    if (--i < 0)
                    {
                        // throw NoSpace("")
                        return false;
                    }

                    Position3 const ipv(draw_ipv(r12, dt(), D12));

                    newpos1 = world_.apply_boundary(
                        particle.position() + ipv * (D1 / D12));
                    newpos2 = world_.apply_boundary(
                        particle.position() - ipv * (D2 / D12));
                    std::vector<
                        std::pair<std::pair<ParticleID, Particle>, Real> >
                        overlapped1(world_.get_particles_within_radius(
                                        newpos1, radius1, pid));
                    std::vector<
                        std::pair<std::pair<ParticleID, Particle>, Real> >
                        overlapped2(world_.get_particles_within_radius(
                                       newpos2, radius2, pid));
                    if (overlapped1.size() == 0 && overlapped2.size() == 0)
                    {
                        break;
                    }
                }

                remove_particle(pid);

                Particle particle_to_update1(
                    species_new1, newpos1, radius1, D1);
                Particle particle_to_update2(
                    species_new2, newpos2, radius2, D2);
                // world_.update_particle(pid1, particle_to_update1);
                // world_.update_particle(pid2, particle_to_update2);
                break;
            }
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
    Real const r12(particle1.radius() + particle2.radius());
    Real const rnd(rng().uniform(0, 1));
    Real prob(0);

    for (ReactionRuleVector::const_iterator i(reaction_rules.begin());
         i != reaction_rules.end(); ++i)
    {
        ReactionRule const& rr(*i);
        prob += rr.k() * dt() / (
            (Igbd_3d(r12, dt(), D1) + Igbd_3d(r12, dt(), D2)) * 4 * M_PI);

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
                remove_particle(pid1);
                remove_particle(pid2);
                break;
            case 1:
                remove_particle(pid1);
                remove_particle(pid2);
                break;
            default:
                throw NotImplemented(
                    "more than one product is not allowed");
                break;
            }
            return true;
        }
    }

    return false;
}

bool BDPropagator::remove_particle(ParticleID const& pid)
{
    world_.remove_particle(pid);
    particle_finder cmp(pid);
    std::vector<std::pair<ParticleID, Particle> >::iterator
        i(std::find_if(queue_.begin(), queue_.end(), cmp));
    if (i != queue_.end())
    {
        queue_.erase(i);
    }
}

} // bd

} // ecell4
