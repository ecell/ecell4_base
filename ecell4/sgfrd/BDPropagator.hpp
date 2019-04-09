#ifndef ECELL4_SGFRD_BD_PROPAGATOR
#define ECELL4_SGFRD_BD_PROPAGATOR
#include <ecell4/core/exceptions.hpp>
#include <ecell4/core/RandomNumberGenerator.hpp>
#include <ecell4/core/Model.hpp>
#include <ecell4/core/geometry.hpp>
#include <ecell4/sgfrd/tracer.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/iterator/iterator_traits.hpp>
#include <boost/typeof/typeof.hpp>
#include <boost/foreach.hpp>
#include <boost/static_assert.hpp>
#include <ecell4/sgfrd/ReactionInfo.hpp>
#include <ecell4/sgfrd/SGFRDWorld.hpp>

namespace ecell4
{

namespace sgfrd
{

/* @brief execute BD algorithm for Multi shell.                      *
 * @tparam containerT is expected to be a World or its wrapper class */
template<typename containerT, typename volume_clearerT>
class BDPropagator
{
public:
    typedef containerT container_type;
    typedef volume_clearerT volume_clearer_type;
    typedef ecell4::Polygon  polygon_type;
    typedef polygon_type::FaceID   FaceID;
    typedef polygon_type::EdgeID   EdgeID;
    typedef polygon_type::VertexID VertexID;

    typedef ecell4::Model   model_type;
    typedef ecell4::Species species_type;
    typedef typename container_type::particle_container_type queue_type;

    // reaction stuff
    typedef ecell4::ReactionRule reaction_rule_type;
    typedef ecell4::sgfrd::ReactionInfo reaction_info_type;
    typedef std::pair<reaction_rule_type, reaction_info_type> reaction_log_type;
    typedef std::vector<reaction_log_type>                reaction_archive_type;

public:

    BDPropagator(const model_type& model, container_type& container,
                 const polygon_type& p, RandomNumberGenerator& rng,
                 const Real dt, const Real rl,
                 reaction_archive_type& last_reactions,
                 volume_clearer_type vc)
    : max_retry_count_(1), dt_(dt), reaction_length_(rl), rejected_move_count_(0),
      container_(container), model_(model), polygon_(p), rng_(rng),
      last_reactions_(last_reactions), vc_(vc), queue_(container.list_particles())
    {
        shuffle(rng_, queue_);
    }

    bool operator()()
    {
        SGFRD_SCOPE(ns, BDPropagator, this->vc_.access_tracer())
        if(queue_.empty())
        {
            return false;
        }

        // make copy of the next particle
        ParticleID pid; Particle p;
        boost::tie(pid, p) = queue_.back(); queue_.pop_back();
        FaceID fid = this->container_.get_face_id(pid);

        // to restore the position, copy previous state.
        const Real3  prev_pos(p.position());
        const FaceID prev_fid(fid);

        if(this->attempt_reaction(pid, p, fid))
        {
            return true;
        }
        if(p.D() == 0.0)
        {
            return true;
        }

        // no 1st order reaction occured & particle is movable.
        auto position     = std::make_pair(p.position(), fid);
        auto displacement = this->draw_displacement(p, fid);
        this->propagate(position, displacement);

        SGFRD_TRACE(this->vc_.access_tracer().write(
                    "particle %1% propagated", pid));

        // check escapement and clear volume if needed
        {
            // update local copy of particle
            boost::tie(p.position(), fid) = position;
            if(!clear_volume(p, fid, pid))
            {
                // rejected. restore position. previous position does not cause
                // overlap because position and species are kept intact.
                ++(this->rejected_move_count_);
                p.position() = prev_pos;
                fid          = prev_fid;
            }
        }

        // retrieve possible reactants (within r1+r2+reaction_length)
        auto overlapped = this->list_reaction_overlap(pid, p, fid);

        // check core-overlap
        std::pair<ParticleID, Particle> pp; Real d;
        bool core_overlapped = false;
        for(const auto& ppd : overlapped)
        {
            std::tie(pp, d) = ppd;
            if(d < p.radius() + pp.second.radius())
            {
                // core overlap!
                // restore position and re-collect overlapped particles
                p.position() = prev_pos;
                fid          = prev_fid;
                core_overlapped = true;
                break;
            }
        }

        if(core_overlapped)
        {
            overlapped = this->list_reaction_overlap(pid, p, fid);
        }
        else
        {
            this->container_.update_particle(pid, p, fid);
        }

        if(overlapped.empty())
        {
            // no reaction-partner exists. overlaps are already cleared. update.
            return true;
        }

        // attempt 2nd order reaction...
        const bool react = this->attempt_reaction(
                pid, p, fid, overlapped.begin(), overlapped.end());
        if(!react)
        {
            ++(this->rejected_move_count_);
        }
        return true;
    }

    Real                       dt()  const throw() {return dt_;}
    RandomNumberGenerator&     rng()       throw() {return rng_;}
    volume_clearer_type const& vc()  const throw() {return vc_;}
    std::size_t rejected_moves() const throw() {return this->rejected_move_count_;}

  protected:

    bool attempt_reaction(const ParticleID& pid, const Particle& p,
                          const FaceID& fid)
    {
        SGFRD_SCOPE(ns, BD_attempt_single_reaction, this->vc_.access_tracer())

        const auto& rules = this->model_.query_reaction_rules(p.species());
        SGFRD_TRACE(this->vc_.access_tracer().write(
                    "%1% rules found for particle %2%", rules.size(), pid))
        if(rules.empty())
        {
            return false;
        }

        const Real rnd(this->rng_.uniform(0., 1.));
        SGFRD_TRACE(this->vc_.access_tracer().write("drawn probability = %1%", rnd))
        Real prob = 0.;
        for(const auto& rule : rules)
        {
            SGFRD_TRACE(this->vc_.access_tracer().write("k * dt = %1%",
                        rule.k() * dt_))
            if((prob += rule.k() * dt_) <= rnd)
            {
                continue;
            }
            if(prob >= 1.0)
            {
                std::cerr << "reaction prob exceeds 1" << std::endl;
            }

            switch(rule.products().size())
            {
                case 0:
                {
                    SGFRD_TRACE(this->vc_.access_tracer().write(
                        "1->0 reaction occured."))
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

    template<typename Iterator>
    bool attempt_reaction(
        const ParticleID& pid1, const Particle& p1, const FaceID& f1,
        const Iterator first, const Iterator last)
    {
        // Iterator::value_type == pair<pair<ParticleID, Particle>, Real>;
        static_assert(boost::is_same<
            typename boost::iterator_value<Iterator>::type,
            std::pair<std::pair<ParticleID, Particle>, Real> >::value, "");

        SGFRD_SCOPE(ns, BD_attempt_pair_reaction, this->vc_.access_tracer())

        const Real rnd(rng_.uniform(0., 1.));
        Real acc_prob = 0.;

        for(Iterator iter(first); iter != last; ++iter)
        {
            const ParticleID& pid2 = iter->first.first;
            const Particle&   p2   = iter->first.second;

            const auto& rules =
                this->model_.query_reaction_rules(p1.species(), p2.species());
            if(rules.empty())
            {
                // no reaction can occur because there is no rule.
                continue;
            }

            const Real k_tot = this->calc_k_total(rules);
            acc_prob += k_tot * calc_acceptance_coef(p1, p2);

            if(acc_prob <= rnd)
            {
                continue;
            }
            else if(1.0 < acc_prob)
            {
                std::cerr << "WARNING: reaction probability exceeds 1\n";
            }

            const auto& rule = this->determine_reaction_rule_from(rules, k_tot);
            switch(rule.products().size())
            {
                case 0: // 2->0 reaction
                {
                    SGFRD_TRACE(this->vc_.access_tracer().write("particle "
                                "%1% and %2% degradated", pid1, pid2));
                    this->remove_particle(pid1);
                    this->remove_particle(pid2);
                    this->last_reactions_.push_back(std::make_pair(rule,
                         init_reaction_info(pid1, p1, pid2, p2)));
                    return true;
                }
                case 1:
                {
                    const FaceID& f2 = this->container_.get_face_id(pid2);
                    const bool reacted = this->attempt_reaction_2_to_1(
                        pid1, p1, f1, pid2, p2, f2, std::make_pair(rule,
                            init_reaction_info(pid1, p1, pid2, p2)));
                    if(reacted)
                    {
                        SGFRD_TRACE(this->vc_.access_tracer().write("particle "
                                    "%1% and %2% bound", pid1, pid2));
                        return true;
                    }
                    else
                    {
                        return false;
                    }
                }
                default:
                {
                    throw NotSupported("BDPropagator: 2 -> N (N>1) "
                                       "reaction is not allowed");
                }
            }
        }
        return false;
    }

    /*! @brief A -> B case.
     * particle does not move. but its radius may change. if overlap occur after
     * reaction, the reaction will be rejected. if accepted, this function does
     * both update and record reaction. */
    bool attempt_reaction_1_to_1(
            const ParticleID& pid, const Particle& p, const FaceID& fid,
            reaction_log_type rlog)
    {
        SGFRD_SCOPE(ns, BD_attempt_1to1_reaction, this->vc_.access_tracer())
        const species_type species_new =
            this->model_.apply_species_attributes(rlog.first.products().front());
        const Real radius_new = species_new.get_attribute_as<Real>("radius");
        const Real D_new      = species_new.get_attribute_as<Real>("D");

        if(is_overlapping(std::make_pair(p.position(), fid), radius_new, pid))
        {
            SGFRD_TRACE(this->vc_.access_tracer().write(
                "1->1 reaction rejected because of the overlapping"))
            return false;
        }

        Particle particle_new(species_new, p.position(), radius_new, D_new);

        if(!clear_volume(particle_new, fid, pid))
        {
            SGFRD_TRACE(this->vc_.access_tracer().write(
                "1->1 reaction rejected because of the overlapping(vc)"))
            return false;
        }

        SGFRD_TRACE(this->vc_.access_tracer().write("1->1 reaction occured"))
        this->container_.update_particle(pid, particle_new, fid);
        rlog.second.add_product(std::make_pair(pid, particle_new));
        last_reactions_.push_back(rlog);
        return true;
    }

    /*! @brief A -> B + C case.
     * after reaction, the products will try to move. */
    bool attempt_reaction_1_to_2(
            const ParticleID& pid, const Particle& p, const FaceID& fid,
            reaction_log_type rlog)
    {
        SGFRD_SCOPE(ns, BD_attempt_1to2_reaction, this->vc_.access_tracer())

        const Species sp1 =
            model_.apply_species_attributes(rlog.first.products().at(0));
        const Species sp2 =
            model_.apply_species_attributes(rlog.first.products().at(1));

        const Real D1  = sp1.get_attribute_as<Real>("D");
        const Real D2  = sp2.get_attribute_as<Real>("D");
        const Real D12 = D1 + D2;
        const Real r1  = sp1.get_attribute_as<Real>("radius");
        const Real r2  = sp2.get_attribute_as<Real>("radius");
        const Real r12 = r1 + r2;

        if(D1 == 0. && D2 == 0)
        {
            throw NotSupported("BDPropagator::1->2: "
                    "reaction between immobile particles");
        }

        const Real3 n = polygon_.triangle_at(fid).normal();

        boost::array<std::pair<Real3, FaceID>, 2> newpfs;
        newpfs[0] = std::make_pair(p.position(), fid);
        newpfs[1] = std::make_pair(p.position(), fid);

        const Real  separation_factor = r12 * 1e-7;
        std::size_t separation_count  = 10u;
        while(separation_count != 0)
        {
            --separation_count;
            SGFRD_TRACE(this->vc_.access_tracer().write(
                        "separation count = %1%", separation_count));

            const Real3 ipv(draw_ipv(r12 + separation_factor, D12, n));
            Real3 disp1(ipv * (D1 / D12)), disp2(ipv * (-D2 / D12));

            // put two particles next to each other
            this->propagate(newpfs[0], disp1);
            this->propagate(newpfs[1], disp2);

            const Real dist = ecell4::polygon::distance(this->polygon_,
                    newpfs[0], newpfs[1]);
            if(dist <= r12) // check the new positions
            {
                newpfs[0] = std::make_pair(p.position(), fid); //rollback
                newpfs[1] = std::make_pair(p.position(), fid);
                continue;
            }

            if(is_overlapping(newpfs[0], r1, pid) ||
               is_overlapping(newpfs[1], r2, pid))
            {
                SGFRD_TRACE(this->vc_.access_tracer().write(
                    "1->2 reaction rejected because of no space"))
                return false; // no space
            }
        }

        boost::array<Particle, 2> particles_new;
        particles_new[0] = Particle(sp1, newpfs[0].first, r1, D1);
        particles_new[1] = Particle(sp2, newpfs[1].first, r2, D2);

        if(!clear_volume(particles_new[0], newpfs[0].second, pid))
        {
            SGFRD_TRACE(this->vc_.access_tracer().write(
                "1->2 reaction rejected because clear_volume failed (no space)"))
            return false;
        }
        if(!clear_volume(particles_new[1], newpfs[1].second, pid))
        {
            SGFRD_TRACE(this->vc_.access_tracer().write(
                "1->2 reaction rejected because clear_volume failed (no space)"))
            return false;
        }

        const bool update_result = this->container_.update_particle(
                pid, particles_new[0], newpfs[0].second);
        auto pp2 =
            this->container_.new_particle(particles_new[1], newpfs[1].second);
        const ParticleID pid2(pp2.first.first);

        assert(update_result == false); // no particle generation occured
        assert(pp2.second    == true);  // new particle is generated

        SGFRD_TRACE(this->vc_.access_tracer().write("1->2 reaction occured"))

        //----------------------------- trial move -----------------------------

        Integer num_move_particle = rng_.uniform_int(1, 2);
        bool move_first_particle = (rng_.uniform_int(0, 1) == 0);

        while(num_move_particle != 0)
        {
            const ParticleID pid_to_move = (move_first_particle) ? pid : pid2;
            const FaceID fid_to_move =
                this->container_.get_face_id(pid_to_move);
            Particle p_to_move =
                this->container_.get_particle(pid_to_move).second;

            std::pair<Real3, FaceID> position = std::make_pair(
                    p_to_move.position(), fid_to_move);
            Real3 displacement = draw_displacement(p_to_move, fid_to_move);
            this->propagate(position, displacement);

            if(!is_overlapping(position, p_to_move.radius(), pid_to_move))
            {
                const Real3 backup = p_to_move.position();
                p_to_move.position() = position.first;
                if(clear_volume(p_to_move, position.second, pid_to_move))
                {
                    this->container_.update_particle(
                            pid_to_move, p_to_move, position.second);
                }
                else
                {
                    p_to_move.position() = backup;
                }
            }

            --num_move_particle;
            move_first_particle = !move_first_particle;
        }

        rlog.second.add_product(
            std::make_pair(pid,  this->container_.get_particle(pid).second));
        rlog.second.add_product(
            std::make_pair(pid2, this->container_.get_particle(pid).second));
        last_reactions_.push_back(rlog);

        return true;
    }

    bool attempt_reaction_2_to_1(// XXX consider using boost::tuple
            const ParticleID& pid1, const Particle& p1, const FaceID& fid1,
            const ParticleID& pid2, const Particle& p2, const FaceID& fid2,
            reaction_log_type rlog)
    {
        const species_type sp_new(rlog.first.products().front());
        const Real radius_new = sp_new.get_attribute_as<Real>("radius");
        const Real D_new      = sp_new.get_attribute_as<Real>("D");

        const Real3 pos1(p1.position()), pos2(p2.position());
        const Real D1(p1.D()), D2(p2.D());
        const Real D12(D1 + D2);

        // this calculates the center position weighted by Diffusion coef.

        std::pair<Real3, FaceID> pf1;
        if(D1 == 0.0)
        {
            pf1 = std::make_pair(pos1, fid1);
        }
        else if(D2 == 0.0)
        {
            pf1 = std::make_pair(pos2, fid2);
        }
        else
        {
            Real3 dp = ecell4::polygon::direction(this->polygon_,
                std::make_pair(pos1, fid1), std::make_pair(pos2, fid2)) * (D1 / D12);
            pf1 = std::make_pair(pos1, fid1);
            this->propagate(pf1, dp);
        }

        if(is_overlapping(pf1, radius_new, pid1, pid2))
        {
            std::cout << "ERROR: the new particle (A+B->B) overlaps!" << std::endl;
            assert(false);
            return false;
        }

        const Particle particle_new(sp_new, pf1.first, radius_new, D_new);
        if(!clear_volume(particle_new, pf1.second, pid1, pid2))
        {
            std::cout << "ERROR: the new particle (A+B->B) overlaps!" << std::endl;
            assert(false);
            return false;
        }

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

    void remove_particle(const ParticleID& pid, const FaceID& fid)
    {
        container_.remove_particle(pid, fid);
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
            const std::pair<Real3, FaceID>& pos, const Real& rad,
            const ParticleID& pid) const
    {
        return !(this->container_.check_no_overlap(pos, rad, pid));
    }
    bool is_overlapping(
            const std::pair<Real3, FaceID>& pos, const Real& rad,
            const ParticleID& pid1, const ParticleID& pid2) const
    {
        return !(this->container_.check_no_overlap(pos, rad, pid1, pid2));
    }

    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
    list_reaction_overlap(const ParticleID& pid, const Particle& p,
                          const FaceID& fid) const
    {
        return this->container_.list_particles_within_radius(
            std::make_pair(p.position(), fid), p.radius() + reaction_length_, pid);
    }

    void propagate(std::pair<Real3, FaceID>& pos, Real3& disp) const
    {
        pos = ecell4::polygon::travel(this->polygon_, pos, disp);
        return ;
    }

    /*! clear the region that particle will occupy. returns if succeed.
     *  this may burst some domains that overlap with particle. */
    bool clear_volume(const Particle& p, const FaceID& fid)
    {
        return vc_(p, fid);
    }

    bool clear_volume(const Particle& p, const FaceID& fid,
                      const ParticleID& ignore)
    {
        return vc_(p, fid, ignore);
    }

    bool clear_volume(const Particle& p, const FaceID& fid,
                      const ParticleID& ignore1, const ParticleID& ignore2)
    {
        return vc_(p, fid, ignore1, ignore2);
    }

    Real3 random_circular_uniform(const Real& r)
    {
        const Real theta = this->rng_.uniform(0., 2 * M_PI);
        return Real3(r * std::cos(theta), r * std::sin(theta), 0.);
    }

    Real3 random_circular_uniform(const Real r, const Real3& normal)
    {
        const Real3 rnd = random_circular_uniform(r);
        const Real tilt = calc_angle(Real3(0, 0, 1), normal);

        if     (std::abs(tilt)        < 1e-10) {return rnd;}
        else if(std::abs(tilt - M_PI) < 1e-10) {return rnd * (-1.0);}
        else if(std::abs(tilt + M_PI) < 1e-10) {return rnd * (-1.0);}

        const Real3 axis = cross_product(Real3(0., 0., 1.), normal);
        return rotate(tilt, axis * (1. / length(axis)), rnd);
    }

    Real3 draw_displacement(const Particle& p, const FaceID& fid)
    {
        const Real r = rng_.gaussian(std::sqrt(4 * p.D() * dt_));
        return random_circular_uniform(r, polygon_.triangle_at(fid).normal());
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
        if((p1.D() == 0) || (p2.D() == 0))
        {
            // immovable particles immediately return after attempting 1st order
            // reaction.
            // to attempt 2nd order reaction with them, we need to double the
            // acceptance coefficient.
            return dt_ / reaction_area;
        }
        else
        {
            // movable particles checks 2nd order reaction. If the both reactants
            // are movable, both particle attempts reaction. So here it halves
            // the acceptance coefficient to avoid double-counting.
            return 0.5 * dt_ / reaction_area;
        }
    }

    Real calc_k_total(const std::vector<reaction_rule_type>& rules) const noexcept
    {
        Real k_tot(0.0);
        if(rules.empty())
        {
            return k_tot;
        }
        if(rules.size() == 1)
        {
            return rules.front().k();
        }
        for(const auto& rule : rules)
        {
            k_tot += rule.k();
        }
        return k_tot;
    }

    reaction_rule_type const&
    determine_reaction_rule_from(const std::vector<reaction_rule_type>& rules,
                                 const Real k_tot) const noexcept
    {
        assert(!rules.empty());
        if(rules.size() == 1)
        {
            return rules.front();
        }
        const Real rnd(rng_.uniform(0.0, 1.0) * k_tot);
        Real k_cumm(0.0);
        for(const auto& rule : rules)
        {
            k_cumm += rule.k();
            if(rnd < k_cumm)
            {
                return rule;
            }
        }
        return rules.back();
    }

    reaction_info_type init_reaction_info(const ParticleID& pid, const Particle& p)
    {
        return reaction_info_type(container_.t(),
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
        return reaction_info_type(container_.t(), reactants, products);
    }


protected:

    Integer max_retry_count_;
    Real    dt_;
    Real    reaction_length_;
    std::size_t rejected_move_count_;
    container_type&        container_;
    model_type const&      model_;
    polygon_type const&    polygon_;
    RandomNumberGenerator& rng_;
    reaction_archive_type& last_reactions_;
    volume_clearer_type    vc_;
    queue_type queue_;
};

} // sgfrd
} // ecell4
#endif /* ECELL4_SGFRD_BD_PROPAGATOR */
