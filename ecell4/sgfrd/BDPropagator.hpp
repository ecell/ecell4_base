#ifndef ECELL4_SGFRD_BD_PROPAGATOR
#define ECELL4_SGFRD_BD_PROPAGATOR
#include <ecell4/core/exceptions.hpp>
#include <ecell4/core/RandomNumberGenerator.hpp>
#include <ecell4/core/Model.hpp>
#include <ecell4/core/geometry.hpp>
#include <boost/foreach.hpp>
#include "Informations.hpp"
#include "SGFRDWorld.hpp"

namespace ecell4
{

namespace sgfrd
{

/* @brief execute BD algorithm for Multi shell.                      *
 * @tparam containerT is expected to be a World or its wrapper class */
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

    typedef ecell4::Model   model_type;
    typedef ecell4::Species species_type;
    typedef typename container_type::particle_container_type queue_type;

    // reaction stuff
    typedef ecell4::ReactionRule reaction_rule_type;
    typedef ecell4::sgfrd::MoleculeInfo molecule_info_type;
    typedef ecell4::sgfrd::ReactionInfo reaction_info_type;
    typedef std::pair<reaction_rule_type, reaction_info_type> reaction_log_type;
    typedef std::vector<reaction_log_type>                reaction_archive_type;

public:

    BDPropagator(const model_type& model, container_type& container,
                 const polygon_type& p, RandomNumberGenerator& rng,
                 const Real dt, const Real rl,
                 reaction_archive_type& last_reactions)
    : max_retry_count_(1), dt_(dt), reaction_length_(rl), model_(model),
      container_(container), polygon_(p), rng_(rng),
      last_reactions_(last_reactions)
    {
        queue_ = container.list_particles();
        shuffle(rng_, queue_);
    }

    bool operator()()
    {
        if(queue_.empty()) return false;

        ParticleID pid; Particle p;
        boost::tie(pid, p) = queue_.back(); queue_.pop_back();
        BOOST_AUTO(fid, this->container_.get_face_id(pid));

        if(this->attempt_reaction(pid, p, fid)) return true;
        if(p.D() == 0.0)                        return true;

        BOOST_AUTO(position,     std::make_pair(p.position(), fid));
        BOOST_AUTO(displacement, draw_displacement(p, fid));

        this->propagate(position, displacement);

        if(is_overlapping(position, p.radius(), pid))
            position = std::make_pair(p.position(), fid); // restore

        boost::tie(p.position(), fid) = position;

        BOOST_AUTO(overlapped, list_reaction_overlap(pid, p, fid));
        switch(overlapped.size())
        {
            case 0: break;
            case 1:
            {
                ParticleID pid2; Particle p2;
                boost::tie(pid2, p2) = overlapped.front().first;
                attempt_reaction(pid, p, fid, pid2, p2, container_.get_face_id(pid2));
                return true; // check and update is done in attempt_reaction().
            }
            default: return true; // reject if more than 2 particles overlap
        }

        if(clear_volume(p, fid, pid)) // return true if volume is cleared
            this->container_.update_particle(pid, p, fid);

        return true;
    }

    Real                   dt() const {return dt_;}
    RandomNumberGenerator& rng()      {return rng_;}

  protected:

    bool attempt_reaction(const ParticleID& pid, const Particle& p,
                          const face_id_type& fid)
    {
        BOOST_AUTO(const& rules, this->model_.query_reaction_rules(p.species()));
        if(rules.empty()) return false;

        const Real rnd(this->rng_.uniform(0., 1.));
        Real prob = 0.;
        BOOST_FOREACH(reaction_rule_type const& rule, rules)
        {
            if((prob += rule.k() * dt_) <= rnd) continue;
            if(prob >= 1.) std::cerr << "reaction prob exceeds 1" << std::endl;

            switch(rule.products().size())
            {
                case 0:
                {
                    remove_particle(pid);
                    last_reactions_.push_back(std::make_pair(
                                rule, init_reaction_info(pid, p)));
                    return true;
                }
                case 1:
                {
                    return attempt_reaction_1_to_1(pid, p, fid, std::make_pair(
                                rule, init_reaction_info(pid, p)));
                }
                case 2:
                {
                    return attempt_reaction_1_to_2(pid, p, fid, std::make_pair(
                                rule, init_reaction_info(pid, p)));
                }
                default: throw NotImplemented("BDPropagator: "
                    "more than two products from one reactant are not allowed");
            }
        }
        return false;
    }

    bool attempt_reaction(
        const ParticleID& pid1, const Particle& p1, const face_id_type& f1,
        const ParticleID& pid2, const Particle& p2, const face_id_type& f2)
    {
        BOOST_AUTO(const& rules,
            this->model_.query_reaction_rules(p1.species(), p2.species()));
        if(rules.empty()) return 0;

        const Real acceptance_coef = calc_acceptance_coef(p1, p2);
        const Real rnd(rng_.uniform(0., 1.));
        Real prob = 0.;
        BOOST_FOREACH(reaction_rule_type const& rule, rules)
        {
            if((prob += rule.k() * acceptance_coef) <= rnd) continue;
            if(prob >= 1.) std::cerr << "reaction prob exceeds 1" << std::endl;

            switch(rule.products().size())
            {
                case 0:
                {
                    remove_particle(pid1);
                    remove_particle(pid2);
                    last_reactions_.push_back(std::make_pair(rule,
                         init_reaction_info(pid1, p1, pid2, p2)));
                    return true;
                }
                case 1:
                {
                    return attempt_reaction_2_to_1(pid1, p1, f1, pid2, p2, f2,
                        std::make_pair(rule, init_reaction_info(pid1, p1, pid2, p2)));
                }
                default: throw NotImplemented("BDPropagator: "
                    "more than one products from two reactants is not allowed");
            }
        }
        return false;
    }

    /*! @brief A -> B case.
     * particle does not move. but its radius may change. if overlap occur after
     * reaction, the reaction rejected. if accepted, this function does both
     * update and record reaction. */
    bool attempt_reaction_1_to_1(
            const ParticleID& pid, const Particle& p, const face_id_type& fid,
            reaction_log_type rlog)
    {
        const species_type species_new =
            this->model_.apply_species_attributes(rlog.first.products().front());
        const molecule_info_type mol_info =
            this->container_.get_molecule_info(species_new);
        const Real radius_new = mol_info.radius;
        const Real D_new      = mol_info.D;

        if(is_overlapping(std::make_pair(p.position(), fid), radius_new, pid))
            return false;

        Particle particle_new(species_new, p.position(), radius_new, D_new);

        if(!clear_volume(particle_new, fid, pid)) return false; // if no space

        this->container_.update_particle(pid, particle_new, fid);
        rlog.second.add_product(std::make_pair(pid, particle_new));
        last_reactions_.push_back(rlog);
        return true;
    }

    /*! @brief A -> B + C case.
     * after reaction, one of the products will try to move. */
    bool attempt_reaction_1_to_2(
            const ParticleID& pid, const Particle& p, const face_id_type& fid,
            reaction_log_type rlog)
    {
        const Species sp1 =
            model_.apply_species_attributes(rlog.first.products().at(0));
        const Species sp2 =
            model_.apply_species_attributes(rlog.first.products().at(1));

        const molecule_info_type mol1 = container_.get_molecule_info(sp1);
        const molecule_info_type mol2 = container_.get_molecule_info(sp2);

        const Real D1(mol1.D),      D2(mol2.D),      D12(mol1.D + mol2.D);
        const Real r1(mol1.radius), r2(mol2.radius), r12(mol1.radius + mol2.radius);
        const Real3 n = polygon_.triangle_at(fid).normal();

        boost::array<std::pair<Real3, face_id_type>, 2> newpfs;
        newpfs[0] = std::make_pair(p.position(), fid);
        newpfs[1] = std::make_pair(p.position(), fid);

        Integer retry(max_retry_count_);
        while(retry > 0)
        {
            --retry;

            const Real3 ipv(draw_ipv(r12, D12, n));
            Real3 disp1(ipv * ( D1 / D12)), disp2(ipv * (-D2 / D12));

            this->propagate(newpfs[0], disp1);
            this->propagate(newpfs[1], disp2);

            if(is_overlapping(newpfs[0], r1, pid)) continue;
            if(is_overlapping(newpfs[1], r2, pid)) continue;

            break; // no overlap!
        }
        if(retry == 0) return false; // no space

        boost::array<Particle, 2> particles_new;
        particles_new[0] = Particle(sp1, newpfs[0].first, r1, D1);
        particles_new[1] = Particle(sp2, newpfs[1].first, r2, D2);

        if(!clear_volume(particles_new[0], newpfs[0].second, pid)) return false;
        if(!clear_volume(particles_new[1], newpfs[1].second     )) return false;

        this->container_.update_particle(pid, particles_new[0], newpfs[0].second);
        BOOST_AUTO(pp2,
            this->container_.new_particle(particles_new[1], newpfs[1].second));
        const ParticleID pid2(pp2.first.first);

        //----------------------------- tryal move -----------------------------

        if(D1 == 0. && D2 == 0)
            throw std::invalid_argument("reaction between immobile particles");

        const std::size_t idx = // select particle that is moved
            (D2 == 0 || (D1 != 0 && rng_.uniform_int(0, 1) == 0)) ? 0 : 1;
        const ParticleID pid_to_move = (idx==0) ? pid : pid2;

        BOOST_AUTO(position,     newpfs[idx]);
        BOOST_AUTO(displacement, draw_displacement(particles_new[idx],
                                                   newpfs[idx].second));
        this->propagate(position, displacement);
        if(!is_overlapping(position, particles_new[idx].radius(), pid_to_move))
        {
            const Real3 tmp = particles_new[idx].position();
            particles_new[idx].position() = position.first;
            if(clear_volume(particles_new[idx], position.second, pid_to_move))
                this->container_.update_particle(
                        pid_to_move, particles_new[idx], position.second);
            else
                particles_new[idx].position() = tmp;
        }

        rlog.second.add_product(std::make_pair(pid,  particles_new[0]));
        rlog.second.add_product(std::make_pair(pid2, particles_new[1]));
        last_reactions_.push_back(rlog);

        return true;
    }

    bool attempt_reaction_2_to_1(// XXX consider using boost::tuple
            const ParticleID& pid1, const Particle& p1, const face_id_type& fid1,
            const ParticleID& pid2, const Particle& p2, const face_id_type& fid2,
            reaction_log_type rlog)
    {
        const species_type     sp_new(rlog.first.products().front());
        const molecule_info_type info(container_.get_molecule_info(sp_new));
        const Real radius_new(info.radius), D_new(info.D);

        const Real3 pos1(p1.position()), pos2(p2.position());
        const Real D1(p1.D()), D2(p2.D());
        const Real D12(D1 + D2);

        // this calculates the center position weighted by Diffusion coef.
        Real3 dp = polygon_.developed_direction(std::make_pair(pos1, fid1),
                std::make_pair(pos2, fid2)) * (D1 / D12);
        std::pair<Real3, face_id_type> pf1(pos1, fid1);
        this->propagate(pf1, dp);

        if(is_overlapping(pf1, radius_new, pid1, pid2)) return false;

        const Particle particle_new(sp_new, pf1.first, radius_new, D_new);
        if(!clear_volume(particle_new, pf1.second, pid1, pid2)) return false;

        remove_particle(pid2);
        remove_particle(pid1);
        const std::pair<std::pair<ParticleID, Particle>, bool> pp_new =
            this->container_.new_particle(particle_new, pf1.second);

        rlog.second.add_product(pp_new.first);
        last_reactions_.push_back(rlog);

        return true;
    }

    void remove_particle(const ParticleID& pid)
    {
        container_.remove_particle(pid);
        const typename std::vector<std::pair<ParticleID, Particle> >::iterator i(
            std::find_if(queue_.begin(), queue_.end(),
                         ecell4::utils::pair_first_element_unary_predicator<
                             ParticleID, Particle>(pid)));
        if(i != queue_.end())
        {
            queue_.erase(i);
        }
        return;
    }

    bool is_overlapping(
            const std::pair<Real3, face_id_type>& pos, const Real& rad,
            const ParticleID& pid) const
    {
        return !(this->container_.check_no_overlap(pos, rad, pid));
    }
    bool is_overlapping(
            const std::pair<Real3, face_id_type>& pos, const Real& rad,
            const ParticleID& pid1, const ParticleID& pid2) const
    {
        return !(this->container_.check_no_overlap(pos, rad, pid1, pid2));
    }

    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
    list_reaction_overlap(const ParticleID& pid, const Particle& p,
                          const face_id_type& fid) const
    {
        return this->container_.list_particles_within_radius(
            std::make_pair(p.position(), fid), p.radius() + reaction_length_, pid);
    }

    void propagate(std::pair<Real3, face_id_type>& pos, Real3& disp) const
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

    /*! clear the region that particle will occupy. returns if succeed.
     *  this may burst some domains that overlap with particle. */
    bool clear_volume(const Particle& p, const face_id_type& fid)
    {// TODO
        return true;
    }

    bool clear_volume(const Particle& p, const face_id_type& fid,
                      const ParticleID& ignore)
    {// TODO
        return true;
    }

    bool clear_volume(const Particle& p, const face_id_type& fid,
                      const ParticleID& ignore1, const ParticleID& ignore2)
    {// TODO
        return true;
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

             if(std::abs(tilt - M_PI) < 1e-10) return rnd;
        else if(std::abs(tilt + M_PI) < 1e-10) return rnd * (-1.0);

        const Real3 axis = cross_product(Real3(0., 0., 1.), normal);
        return rotate(tilt, axis * (1. / length(axis)), rnd);
    }

    Real3 draw_displacement(const Particle& p, const face_id_type& fid)
    {
        return random_circular_uniform(rng_.gaussian(std::sqrt(4*p.D()*dt_)),
                                       polygon_.triangle_at(fid).normal());
    }

    Real3 draw_ipv(const Real r, const Real D, const Real3& normal)
    {
        const Real rl    = r + this->reaction_length_;
        const Real r_sq  = r * r;
        const Real rd    = rl * rl - r_sq;
        const Real ipvl  = std::sqrt(r_sq + this->rng_.uniform(0., 1.) * rd);
        return random_circular_uniform(ipvl, normal);
    }

    Real calc_reaction_area(const Real radius_sum) const
    {
        const Real rad_react(radius_sum + this->reaction_length_);
        return M_PI * (rad_react * rad_react - radius_sum * radius_sum);
    }

    Real calc_acceptance_coef(const Particle& p1, const Particle& p2) const
    {
        const Real reaction_area = calc_reaction_area(p1.radius() + p2.radius());
        const bool double_count  = ((p1.D()==0) || (p2.D()==0));
        return double_count ? (dt_ / reaction_area) : (0.5 * dt_ / reaction_area);
    }

    reaction_info_type init_reaction_info(const ParticleID& pid, const Particle& p)
    {
        return reaction_info_type(container_.t() + dt_,
                reaction_info_type::container_type(1, std::make_pair(pid, p)),
                reaction_info_type::container_type());
    }

    reaction_info_type
    init_reaction_info(const ParticleID& pid1, const Particle& p1,
                       const ParticleID& pid2, const Particle& p2)
    {
        typename reaction_info_type::container_type reactants(2);
        typename reaction_info_type::container_type products;
        reactants[0] = std::make_pair(pid1, p1);
        reactants[1] = std::make_pair(pid2, p2);
        return reaction_info_type(container_.t() + dt_, reactants, products);
    }


protected:

    Integer max_retry_count_;
    Real    dt_;
    Real    reaction_length_;
    container_type&        container_;
    model_type const&      model_;
    polygon_type const&    polygon_;
    RandomNumberGenerator& rng_;
    reaction_archive_type& last_reactions_;
    queue_type queue_;
};

} // sgfrd
} // ecell4
#endif /* ECELL4_SGFRD_BD_PROPAGATOR */
