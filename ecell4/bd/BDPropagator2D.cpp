#include <iterator>

#include <ecell4/core/exceptions.hpp>
#include <ecell4/core/Species.hpp>

#include "BDPropagator2D.hpp"


namespace ecell4
{

namespace bd
{

bool BDPropagator2D::operator()()
{
    if (queue_.empty())
    {
        return false;
    }

    const ParticleID pid(queue_.back().first);
    queue_.pop_back();
    Particle particle(world_.get_2D_particle(pid).second);
    BDPolygon::face_id_type fid = world_.face_id_on(pid);
    Triangle const& face = world_.face_on(pid);

    if(attempt_reaction(pid, particle, fid))
    {
        return true;
    }

    const Real D(particle.D());
    if(D == 0)
    {
        return true;
    }

    const std::pair<Real3, BDPolygon::face_id_type> newpos(
        world_.apply_surface(std::make_pair(particle.position(), fid),
                             draw_displacement(particle, face.normal())));

    Particle particle_to_update(
        particle.species(), newpos.first, particle.radius(), particle.D());
    // Particle particle_to_update(
    //     particle.species_serial(), newpos, particle.radius(), particle.D());
    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
        overlapped(world_.list_2D_particles_within_radius(
                       newpos, particle.radius(), pid));

    switch (overlapped.size())
    {
    case 0:
        world_.update_2D_particle_without_checking(pid, particle_to_update, newpos.second);
        return true;
    case 1:
        {
            std::pair<ParticleID, Particle, face_id_type> closest(
                (*(overlapped.begin())).first);
            if (attempt_reaction(
                    pid, particle_to_update, fid, closest.first, closest.second))
            {
                return true;
            }
        }
        return true;
    default:
        return true;
    }
}

bool BDPropagator2D::attempt_reaction(
    const ParticleID& pid, const Particle& particle, const BDPolygon::face_id_type& fid)
{
    std::vector<ReactionRule> const& reaction_rules = 
        model_.query_reaction_rules(particle.species());

    if (reaction_rules.size() == 0)
        return false;

    const Real rnd(rng().uniform(0., 1.));
    Real prob = 0.;
    for(std::vector<ReactionRule>::const_iterator
            iter = reaction_rules.begin(); iter != reaction_rules.end(); ++iter)
    {
        const ReactionRule& rule(*iter);
        prob += rule.k() * dt();
        if(prob <= rnd) continue;

        ReactionRule::product_container_type const& products = rule.products();
        reaction_info_type r_info(world_.t() + dt_,
            reaction_info_type::container_type(1, std::make_pair(pid, particle)),
            reaction_info_type::container_type());

        switch(products.size())
        {
        case 0: // decay reaction 1->0
            remove_particle(pid);
            last_reactions_.push_back(std::make_pair(rule, r_info));
            return true;

        case 1: // transform reaction 1->1
            const Species species_new =
                model_.apply_species_attributes(products.front());
            const BDWorld::molecule_info_type mol_info =
                world_.get_molecule_info(species_new);
            const Real radius_new = mol_info.radius;
            const Real D_new      = mol_info.D;

            // confirm whether does new molecule overlaps
            std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
                overlapped = world_.list_2D_particles_within_radius(
                    particle.position(), radius_new, pid, fid);
            if(overlapped.size() != 0)
                return false;

            Particle particle_to_update(
                    species_new, particle.position(), radius_new, D_new);

            world_.update_2D_particle_without_checking(// already checked
                    std::make_pair(pid, particle_to_update, fid));
            r_info.add_product(std::make_pair(pid, particle_to_update));
            last_reactions_.push_back(std::make_pair(rule, r_info));
            return true;

        case 2: // split reaction 1->2
            // TODO
            // XXX: 2D impl of degradation is different from 3D!
            return false;
            break;

        default:
            throw NotImplemented(
                    "BDPropagator: more than two products are not allowed");
        }
    }
    return false;
}

bool BDPropagator2D::attempt_reaction(
    const ParticleID& pid1, const Particle& particle1,
    const ParticleID& pid2, const Particle& particle2)
{
//     std::vector<ReactionRule> reaction_rules(
//         model_.query_reaction_rules(
//             particle1.species(), particle2.species()));
//     if (reaction_rules.size() == 0)
//     {
//         return false;
//     }
//
//     const Real D1(particle1.D()), D2(particle2.D());
//     const Real r12(particle1.radius() + particle2.radius());
//     const Real rnd(rng().uniform(0, 1));
//     Real prob(0);
//
//     for (std::vector<ReactionRule>::const_iterator i(reaction_rules.begin());
//          i != reaction_rules.end(); ++i)
//     {
//         const ReactionRule& rr(*i);
//         prob += rr.k() * dt() / (
//             (Igbd_3d(r12, dt(), D1) + Igbd_3d(r12, dt(), D2)) * 4 * M_PI);
//
//         if (prob >= 1)
//         {
//             // throw std::runtime_error(
//             //     "the total reaction probability exceeds 1."
//             //     " the step interval is too long");
//             std::cerr <<
//                 "the total reaction probability exceeds 1."
//                 " the step interval is too long" << std::endl;
//         }
//         if (prob > rnd)
//         {
//             const ReactionRule::product_container_type& products(rr.products());
//             reaction_info_type ri(world_.t() + dt_, reaction_info_type::container_type(1, std::make_pair(pid1, particle1)), reaction_info_type::container_type());
//             ri.add_reactant(std::make_pair(pid2, particle2));
//
//             switch (products.size())
//             {
//             case 0:
//                 remove_particle(pid1);
//                 remove_particle(pid2);
//
//                 last_reactions_.push_back(std::make_pair(rr, ri));
//                 break;
//             case 1:
//                 {
//                     const Species sp(*(products.begin()));
//                     const BDWorld::molecule_info_type
//                         info(world_.get_molecule_info(sp));
//                     const Real radius_new(info.radius);
//                     const Real D_new(info.D);
//
//                     const Real3 pos1(particle1.position());
//                     const Real3 pos2(
//                         world_.periodic_transpose(particle2.position(), pos1));
//                     const Real D1(particle1.D()), D2(particle2.D());
//                     const Real D12(D1 + D2);
//                     const Real3 newpos(
//                         world_.apply_boundary((pos1 * D2 + pos2 * D1) / D12));
//
//                     std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
//                         overlapped(world_.list_particles_within_radius(
//                                        newpos, radius_new, pid1, pid2));
//                     if (overlapped.size() > 0)
//                     {
//                         // throw NoSpace("");
//                         return false;
//                     }
//
//                     const Particle particle_to_update(
//                         sp, newpos, radius_new, D_new);
//                     remove_particle(pid2);
//                     // world_.update_particle(pid1, particle_to_update);
//                     remove_particle(pid1);
//                     std::pair<std::pair<ParticleID, Particle>, bool> retval = world_.new_particle(particle_to_update);
//
//                     ri.add_product(retval.first);
//                     last_reactions_.push_back(std::make_pair(rr, ri));
//                 }
//                 break;
//             default:
//                 throw NotImplemented(
//                     "more than one product is not allowed");
//                 break;
//             }
//             return true;
//         }
//     }
//
    return false;
}

void BDPropagator2D::remove_particle(const ParticleID& pid)
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
