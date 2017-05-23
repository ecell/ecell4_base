#ifndef ECELL4_SGFRD_BD_PROPAGATOR
#define ECELL4_SGFRD_BD_PROPAGATOR
#include <ecell4/core/RandomNumberGenerator.hpp>
#include <ecell4/core/Model.hpp>
#include <ecell4/core/geometry.hpp>
#include "Informations.hpp"
#include "SGFRDWorld.hpp"

namespace ecell4
{

namespace sgfrd
{

template<typename containerT>
class BDPropagator
{
public:
    typedef containerT container_type;
    typedef ecell4::sgfrd::polygon_traits   polygon_traits_type;
    typedef Polygon<polygon_traits_type>    polygon_type;
    typedef polygon_type::triangle_type     triangle_type;
    typedef polygon_type::face_id_type      face_id_type;
    typedef polygon_type::edge_id_type      edge_id_type;
    typedef polygon_type::vertex_id_type    vertex_id_type;
    typedef polygon_type::face_descripter   face_descripter;
    typedef polygon_type::edge_descripter   edge_descripter;
    typedef polygon_type::vertex_descripter vertex_descripter;
    typedef polygon_type::local_index_type  local_index_type;
    typedef polygon_type::barycentric_type  barycentric_type;

    typedef ecell4::Model model_type;
    typedef typename container_type::particle_container_type queue_type;

    // reaction stuff
    typedef ecell4::ReactionRule reaction_rule_type;
    typedef ecell4::sgfrd::MoleculeInfo molecule_info_type;
    typedef ecell4::sgfrd::ReactionInfo reaction_info_type;
    typedef std::vector<std::pair<ReactionRule, reaction_info_type>
            > reaction_archive_type;

public:

    BDPropagator(Model& model, container_type& container, const polygon_type& p,
                 RandomNumberGenerator& rng, const Real dt, const Real rl,
                 reaction_archive_type& last_reactions)
    : max_retry_count_(1), dt_(dt), reaction_length_(rl), model_(model),
      container_(container), polygon_(p), rng_(rng),
      last_reactions_(last_reactions)
    {
        queue_ = container.list_particles();
        shuffle(rng_, queue_);
    }

    bool operator()();

    Real                   dt() const {return dt_;}
    RandomNumberGenerator& rng()      {return rng_;}

    bool attempt_reaction(const ParticleID& pid, const Particle& particle,
                          const face_id_type& fid);
    bool attempt_reaction(
        const ParticleID& pid1, const Particle& particle1, const face_id_type& f1,
        const ParticleID& pid2, const Particle& particle2, const face_id_type& f2);

    void remove_particle(const ParticleID& pid)
    {
        const typename std::vector<std::pair<ParticleID, Particle> >::iterator i(
            std::find_if(queue_.begin(), queue_.end(),
                         ecell4::utils::pair_first_element_unary_predicator<
                             ParticleID, Particle>(pid)));
        if(i != queue_.end())
        {
            queue_.erase(i);
        }
    }

private:

    void propagate_on_surface(std::pair<Real3, face_id_type>& pos, Real3& disp) const
    {
        unsigned int continue_count = 100;
        while(continue_count > 0)
        {
            boost::tie(pos, disp) = polygon_.move_next_face(pos, disp);
            if(disp[0] == 0. && disp[1] == 0. && disp[2] == 0.) break;
            --continue_count;
        }
        if(continue_count == 0)
            std::cerr << "[WARNING] moving on face by BD: precision lost\n";
        return ;
    }

    Real3 random_circular_uniform(const Real& r)
    {
        const Real theta = this->rng_.uniform(0., 2 * M_PI);
        return Real3(r * std::cos(theta), r * std::sin(theta), 0.);
    }

    Real3 random_circular_uniform(const Real r, const Real3& normal)
    {
        const Real3 rnd = random_circular_uniform(r);
        const Real tilt = angle(Real3(0, 0, 1), normal);

        if(std::abs(tilt - M_PI) < 1e-12)
            return rnd;
        else if(std::abs(tilt + M_PI) < 1e-12)
            return rnd * (-1.0);
        else
            return rotate(tilt, cross_product(Real3(0., 0., 1.), normal), rnd);
    }

    Real3 draw_displacement(const Particle& p, const Real3& normal)
    {
        return random_circular_uniform(
                this->rng_.gaussian(std::sqrt(4 * p.D() * dt_)), normal);
    }

    Real3 random_ipv_2d(const Real r, const Real D, const Real3& normal)
    {
        const Real rl    = r + this->reaction_length_;
        const Real r_sq  = r * r;
        const Real rl_sq = rl * rl;
        const Real ipvl  =
            std::sqrt(r_sq + this->rng_.uniform(0., 1.) * (rl_sq - r_sq));
        return random_circular_uniform(ipvl, normal);
    }

    Real3 draw_ipv(const Real r, const Real D, const Real3& normal)
    {
        return random_ipv_2d(r, D, normal);
    }

    Real calc_reaction_area(const Real radius_sum) const
    {
        const Real rad_react(radius_sum + this->reaction_length_);
        return M_PI * (rad_react * rad_react - radius_sum * radius_sum);
    }

protected:

    Integer max_retry_count_;
    Real    dt_;
    Real    reaction_length_;
    ecell4::Model&         model_;
    container_type&        container_;
    polygon_type const&    polygon_;
    RandomNumberGenerator& rng_;
    reaction_archive_type& last_reactions_;
    queue_type queue_;
};

template<typename containerT>
bool BDPropagator<containerT>::operator()()
{
    if(queue_.empty())
    {
        return false;
    }

    ParticleID pid;
    Particle   particle;
    boost::tie(pid, particle) = queue_.back();
    queue_.pop_back();

    const face_id_type   fid  = container_.get_face_id(pid);
    const triangle_type& face = polygon_.triangle_at(fid);

    if(attempt_reaction(pid, particle, fid))
    {
        return true;
    }

    if(particle.D() == 0.0)
    {
        // reaction for immobile particle is considered
        // by acceptance coefficient in attempt_reaction(p1, p2)
        return true;
    }

    // move
    std::pair<Real3, face_id_type> new_position =
        std::make_pair(particle.position(), fid);
    Real3 displacement = this->draw_displacement(particle, face.normal());

    this->propagate_on_surface(new_position, displacement);

    // confirm whether core overlapped or not
    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
        core_overlapped(container_.list_particles_within_radius(
                        new_position, particle.radius(), pid));

    if(!core_overlapped.empty()) // reject
    {
        new_position = std::make_pair(particle.position(), fid);
    }

    particle.position() = new_position.first;

    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
        overlapped(container_.list_particles_within_radius(
            new_position, particle.radius() + this->reaction_length_, pid));

    switch (overlapped.size())
    {
    case 0:
        container_.update_particle(pid, particle, new_position.second);
        return true;
    case 1:
        {
            std::pair<std::pair<ParticleID, Particle>, Real> closest =
                overlapped.front();
            const face_id_type closest_fid =
                container_.get_face_id(closest.first.first);

            if(attempt_reaction(pid, particle, new_position.second,
                closest.first.first, closest.first.second, closest_fid))
            {
                return true;
            }
        }
        return true;
    default:
        return true;
    }
}

template<typename containerT>
bool BDPropagator<containerT>::attempt_reaction(
    const ParticleID& pid, const Particle& particle, const face_id_type& fid)
{
    std::vector<ecell4::ReactionRule> const& reaction_rules =
        model_.query_reaction_rules(particle.species());

    if (reaction_rules.size() == 0)
        return false;

    const Real rnd(rng_.uniform(0., 1.));
    Real prob = 0.;
    for(typename std::vector<ReactionRule>::const_iterator
        iter = reaction_rules.begin(); iter != reaction_rules.end(); ++iter)
    {
        const ReactionRule& rule(*iter);
        prob += rule.k() * dt();
        if(prob <= rnd) continue;

        // reaction occured!
        typename ReactionRule::product_container_type const& products = rule.products();
        reaction_info_type r_info(container_.t() + dt_,
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
                const molecule_info_type mol_info =
                    container_.get_molecule_info(species_new);
                const Real radius_new = mol_info.radius;
                const Real D_new      = mol_info.D;

                // confirm whether does new molecule overlaps
                std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
                    overlapped = container_.list_particles_within_radius(
                        std::make_pair(particle.position(), fid), radius_new,
                        /*ignore1 = */ pid);
                if(overlapped.size() != 0)
                    return false;

                Particle particle_to_update(
                        species_new, particle.position(), radius_new, D_new);
                container_.update_particle(pid, particle_to_update, fid);

                r_info.add_product(std::make_pair(pid, particle_to_update));
                last_reactions_.push_back(std::make_pair(rule, r_info));
                return true;
            }
            case 2: // split reaction 1->2
            {
                const Species sp1 = model_.apply_species_attributes(products.at(0));
                const Species sp2 = model_.apply_species_attributes(products.at(1));
                const molecule_info_type mol1 = container_.get_molecule_info(sp1);
                const molecule_info_type mol2 = container_.get_molecule_info(sp2);

                const Real D1 = mol1.D;
                const Real D2 = mol2.D;
                const Real D12 = D1 + D2;

                const Real r1 = mol1.radius;
                const Real r2 = mol2.radius;
                const Real r12 = r1 + r2;

                const Real3 n = polygon_.triangle_at(fid).normal();

                Integer retry(max_retry_count_);
                boost::array<std::pair<Real3, face_id_type>, 2> newpfs;
                newpfs[0] = std::make_pair(particle.position(), fid);
                newpfs[1] = std::make_pair(particle.position(), fid);

                while(true)
                {
                    if(--retry < 0)
                    {
                        // dissociation rejected because there is no space
                        return false;
                    }
                    const Real3 ipv(draw_ipv(r12, D12, n));
                    Real3 disp1(ipv * ( D1 / D12));
                    Real3 disp2(ipv * (-D2 / D12));

                    this->propagate_on_surface(newpfs[0], disp1);
                    this->propagate_on_surface(newpfs[1], disp2);

                    const std::vector<
                        std::pair<std::pair<ParticleID, Particle>, Real> >
                        overlapped1(container_.list_particles_within_radius(
                            newpfs[0], r1, pid));
                    if(overlapped1.size() > 0) continue;

                    const std::vector<
                        std::pair<std::pair<ParticleID, Particle>, Real> >
                        overlapped2(container_.list_particles_within_radius(
                            newpfs[1], r2, pid));
                    if(overlapped2.size() == 0) break;
                }

                boost::array<Particle, 2> particles_to_update;
                particles_to_update[0] = Particle(sp1, newpfs[0].first, r1, D1);
                particles_to_update[1] = Particle(sp2, newpfs[1].first, r2, D2);

                const bool update_result = container_.update_particle(
                        pid, particles_to_update[0], newpfs[0].second);
                assert(!update_result);

                std::pair<std::pair<ParticleID, Particle>, bool>
                    p2_added = container_.new_particle(
                            particles_to_update[1], newpfs[1].second);
                assert(p2_added.second);
                const ParticleID p2_id(p2_added.first.first);

                // move one particle.
                // if rejected, particles are left at the position
                // just after dissociation reaction.

                assert(D1 != 0 || D2 != 0);
                const std::size_t idx = // select particle that is moved
                    (D2 == 0 || (D1 != 0 && rng_.uniform_int(0, 1) == 0)) ? 0 : 1;

                const Real3& normal =
                    this->polygon_.triangle_at(newpfs[idx].second).normal();

                std::pair<Real3, face_id_type> newpf(newpfs[idx]);
                Real3 disp = draw_displacement(particles_to_update[idx], normal);
                this->propagate_on_surface(newpf, disp);
                std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
                    overlapped(container_.list_particles_within_radius(
                            newpf, particles_to_update[idx].radius(), pid));

                if(overlapped.size() == 0)
                {
                    particles_to_update[idx].position() = newpf.first;
                    const bool update_1 = container_.update_particle(
                            pid, particles_to_update[idx], newpf.second);
                    assert(!update_1);
                }

                r_info.add_product(std::make_pair(pid,   particles_to_update[0]));
                r_info.add_product(std::make_pair(p2_id, particles_to_update[1]));
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

template<typename containerT>
bool BDPropagator<containerT>::attempt_reaction(
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
    const Real reaction_area = calc_reaction_area(r12);
    const Real D1(particle1.D()), D2(particle2.D());
    const Real rnd(rng().uniform(0, 1));
    const bool double_count = particle1.D() == 0 || particle2.D() == 0;
    const Real coef_acceptance_prob = (double_count) ?
        (dt() / reaction_area) : (0.5 * dt() / reaction_area);
    Real prob(0);

    for(typename std::vector<ReactionRule>::const_iterator
            i(reaction_rules.begin()); i != reaction_rules.end(); ++i)
    {
        const ReactionRule& rr(*i);
        prob += rr.k() * coef_acceptance_prob;

        if(prob <= rnd)
        {
            continue;
        }
        else if (prob >= 1)
        {
            std::cerr << "the total reaction probability exceeds 1. "
                      << "the step interval is too long or "
                      << "reaction length is too short" << std::endl;
        }

        // reaction occured (1 > prob > rnd)

        const typename ReactionRule::product_container_type& products(rr.products());
        reaction_info_type ri(container_.t() + dt_,
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
                const molecule_info_type info(container_.get_molecule_info(sp_new));
                const Real radius_new(info.radius);
                const Real D_new(info.D);

                const Real3 pos1(particle1.position());
                const Real3 pos2(particle2.position());
                const Real D1(particle1.D()), D2(particle2.D());
                const Real D12(D1 + D2);

                // this calculates the center position waited by Diffusion coef.
                const Real3 dp = polygon_.developed_direction(
                        std::make_pair(pos1, f1), std::make_pair(pos2, f2));
                Real3 dp1 = dp * (D1 / D12);
                std::pair<Real3, face_id_type> pf1(pos1, f1);
                this->propagate_on_surface(pf1, dp1);

                std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
                    overlapped(container_.list_particles_within_radius(
                                   pf1, radius_new, pid1, pid2));
                if(overlapped.size() > 0)
                {
                    return false;
                }

                const Particle particle_to_update(
                        sp_new, pf1.first, radius_new, D_new);
                remove_particle(pid2);
                remove_particle(pid1);

                std::pair<std::pair<ParticleID, Particle>, bool> retval =
                    container_.new_particle(particle_to_update, pf1.second);

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

} // sgfrd
} // ecell4
#endif /* ECELL4_SGFRD_BD_PROPAGATOR */
