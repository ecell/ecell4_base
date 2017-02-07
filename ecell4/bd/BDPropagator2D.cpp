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

    ParticleContainer2D& container2D = world_.container_2D();

    const ParticleID pid(queue_.back().first);
    queue_.pop_back();
    Particle        particle = container2D.get_particle(pid).second;
    face_id_type    fid = container2D.belonging_faceid(pid);
    Triangle const& face = container2D.belonging_face(pid);

    if(attempt_reaction(pid, particle, fid))
    {
        return true;
    }

    const Real D(particle.D());
    if(D == 0.0)
    {
        return true; //XXX consider reaction
    }

    const std::pair<Real3, face_id_type> newpos(world_.apply_surface(
        std::make_pair(particle.position(), fid),
        draw_displacement(particle, face.normal())));

    Particle particle_to_update(
        particle.species(), newpos.first, particle.radius(), particle.D());

    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
        overlapped(container2D.list_particles_within_radius(
                    newpos, particle.radius() + reaction_length_, pid));

    switch (overlapped.size())
    {
    case 0:
        container2D.update_particle(pid, particle_to_update, newpos.second);
        return true;
    case 1:
        {
            std::pair<std::pair<ParticleID, Particle>, face_id_type> closest =
                overlapped.front();
            if (attempt_reaction(
                    pid, particle_to_update, fid,
                    closest.first.first, closest.first.second, closest.second))
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
    const ParticleID& pid, const Particle& particle, const face_id_type& fid)
{
    std::vector<ReactionRule> const& reaction_rules =
        model_.query_reaction_rules(particle.species());

    if (reaction_rules.size() == 0)
        return false;

    const Real rnd(rng_.uniform(0., 1.));
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
            {
                remove_particle(pid);
                last_reactions_.push_back(std::make_pair(rule, r_info));
                return true;
            }
            case 1: // transform reaction 1->1
            {
                const Species species_new =
                    model_.apply_species_attributes(products.front());
                const BDWorld::molecule_info_type mol_info =
                    world_.get_molecule_info(species_new);
                const Real radius_new = mol_info.radius;
                const Real D_new      = mol_info.D;

                // confirm whether does new molecule overlaps
                std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
                    overlapped = world_.list_particles_within_radius(
                        std::make_pair(particle.position(), fid), radius_new,
                        /*ignore1 = */ pid);
                if(overlapped.size() != 0)
                    return false;

                Particle particle_to_update(
                        species_new, particle.position(), radius_new, D_new);
                world_.update_particle_without_checking(
                        pid, particle_to_update, fid);

                r_info.add_product(std::make_pair(pid, particle_to_update));
                last_reactions_.push_back(std::make_pair(rule, r_info));
                return true;
            }
            case 2: // split reaction 1->2
            {
                const Species sp1 = model_.apply_species_attributes(products.at(0));
                const Species sp2 = model_.apply_species_attributes(products.at(1));
                const molecule_info_type mol1 = world_.get_molecule_info(sp1);
                const molecule_info_type mol2 = world_.get_molecule_info(sp2);

                const Real D1 = mol1.D;
                const Real D2 = mol2.D;
                const Real D12 = D1 + D2;

                const Real r1 = mol1.radius;
                const Real r2 = mol2.radius;
                const Real r12 = r1 + r2;

                const Real3 n = this->poly_.at(fid).normal();

                Integer retry(max_retry_count_);
                std::pair<Real3, face_id_type> newpf1, newpf2;
                while(true)
                {
                    if(--retry < 0)
                    {
                        // dissociation rejected because there is no space
                        return false;
                    }
                    const Real3 ipv(draw_ipv(r12, D12, n));

                    newpf1 = world_.apply_surface(
                            std::make_pair(particle.position(), fid),
                            ipv * (D1 / D12));

                    newpf2 = world_.apply_surface(
                            std::make_pair(particle.position(), fid),
                            ipv * (-D2 / D12));

                    const std::vector<
                        std::pair<std::pair<ParticleID, Particle>, Real> >
                        overlapped1(world_.list_particles_within_radius(
                            newpf1, r1, pid));
                    if(overlapped1.size() != 0) continue;

                    const std::vector<
                        std::pair<std::pair<ParticleID, Particle>, Real> >
                        overlapped2(world_.list_particles_within_radius(
                            newpf2, r2, pid));
                    if(overlapped2.size() == 0) break;
                }

                Particle particle_to_update1(sp1, newpf1.first, r1, D1);
                Particle particle_to_update2(sp2, newpf2.first, r2, D2);

                // move.
                // if rejected, update particles with positions just after reaction

                assert(D1 != 0 || D2 != 0);
                if(D2 == 0 || (D1 != 0 && rng_.uniform_int(0, 1) == 0))
                {
                    const Real3& normal = this->poly_.at(newpf1.second).normal();
                    const ParticleContainer2D& container2D = world_.container_2D();

                    std::pair<Real3, face_id_type> newpf(world_.apply_surface(
                        newpf1, draw_displacement(particle_to_update1, normal)));

                    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
                        overlapped(container2D.list_particles_within_radius(
                                newpf, particle.radius() + reaction_length_, pid));

                    if(overlapped.size() == 0)
                    {
                        particle_to_update1.position() = newpf.first;
                        newpf1.second                  = newpf.second;
                    }
                }
                else
                {// particle 2
                    const Real3& normal = this->poly_.at(newpf2.second).normal();
                    const ParticleContainer2D& container2D = world_.container_2D();

                    std::pair<Real3, face_id_type> newpf(world_.apply_surface(
                        newpf2, draw_displacement(particle_to_update2, normal)));

                    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
                        overlapped(container2D.list_particles_within_radius(
                                newpf, particle.radius() + reaction_length_, pid));

                    if(overlapped.size() == 0)
                    {
                        particle_to_update2.position() = newpf.first;
                        newpf2.second                  = newpf.second;
                    }
                }

                world_.update_particle_without_checking(
                        pid, particle_to_update1, newpf1.second);
                const std::pair<std::pair<ParticleID, Particle>, bool> retval =
                    world_.new_particle(particle_to_update2, newpf2.second);

                r_info.add_product(std::make_pair(pid, particle_to_update1));
                r_info.add_product(retval.first);
                last_reactions_.push_back(std::make_pair(rule, r_info));

                return true;
            }
            default:
            {
                throw NotImplemented(
                        "BDPropagator: more than two products are not allowed");
            }
        }
    }
    return false;
}

bool BDPropagator2D::attempt_reaction(
        const ParticleID& pid1, const Particle& particle1, const face_id_type& f1,
        const ParticleID& pid2, const Particle& particle2, const face_id_type& f2)
{
    std::vector<ReactionRule> reaction_rules(
        model_.query_reaction_rules(
            particle1.species(), particle2.species()));
    if (reaction_rules.size() == 0)
    {
        return false;
    }

    const Real r12(particle1.radius() + particle2.radius());

    const Real D1(particle1.D()), D2(particle2.D());
    const Real rnd(rng().uniform(0, 1));
    Real prob(0);

    for (std::vector<ReactionRule>::const_iterator i(reaction_rules.begin());
         i != reaction_rules.end(); ++i)
    {
        const ReactionRule& rr(*i);
        prob += rr.k() * dt(); //XXX:INCORRECT consider reaction length

        if(prob <= rnd)
        {
            continue;
        }
        else if (prob >= 1)
        {
            std::cerr << "the total reaction probability exceeds 1."
                      << " the step interval is too long" << std::endl;
        }

        // reaction occured (1 > prob > rnd)

        const ReactionRule::product_container_type& products(rr.products());
        reaction_info_type ri(world_.t() + dt_,
                              reaction_info_type::container_type(
                                  1, std::make_pair(pid1, particle1)),
                                  reaction_info_type::container_type());
        ri.add_reactant(std::make_pair(pid2, particle2));

        switch (products.size())
        {
            case 0: // degradation reaction: 2 -> 0
            {
                remove_particle(pid1);
                remove_particle(pid2);

                last_reactions_.push_back(std::make_pair(rr, ri));
                return true;
            }
            case 1: // binding reaction : 2 -> 1
            {
                const Species sp_new(products.front());
                const BDWorld::molecule_info_type
                    info(world_.get_molecule_info(sp_new));
                const Real radius_new(info.radius);
                const Real D_new(info.D);

                const Real3 pos1(particle1.position());
                const Real3 pos2(particle2.position());
                const Real D1(particle1.D()), D2(particle2.D());
                const Real D12(D1 + D2);

                // this calculates the center position waited by Diffusion coef.
                const Real3 dp = world_.get_inter_position_vector(
                        std::make_pair(pos1, f1), std::make_pair(pos2, f2));
                const Real3 dp1 = dp * (D1 / D12);
                const std::pair<Real3, face_id_type> newpf(
                        world_.apply_surface(std::make_pair(pos1, f1), dp1));

                std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
                    overlapped(world_.list_particles_within_radius(
                                   newpf, radius_new, pid1, pid2));
                if(overlapped.size() > 0)
                {
                    return false;
                }

                const Particle particle_to_update(
                        sp_new, newpf.first, radius_new, D_new);
                remove_particle(pid2);
                remove_particle(pid1);

                std::pair<std::pair<ParticleID, Particle>, bool> retval =
                    world_.new_particle(particle_to_update, newpf.second);

                ri.add_product(retval.first);
                last_reactions_.push_back(std::make_pair(rr, ri));
                return true;
            }
            default:
                throw NotImplemented("more than one product is not allowed");
        }
    }
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
