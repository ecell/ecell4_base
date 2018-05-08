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
#include <ecell4/sgfrd/Informations.hpp>
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
    typedef ecell4::Triangle triangle_type;
    typedef polygon_type::FaceID   FaceID;
    typedef polygon_type::EdgeID   EdgeID;
    typedef polygon_type::VertexID VertexID;

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
                 reaction_archive_type& last_reactions,
                 volume_clearer_type vc)
    : max_retry_count_(1), dt_(dt), reaction_length_(rl), rejected_move_count_(0),
      model_(model), container_(container), polygon_(p), rng_(rng),
      last_reactions_(last_reactions), vc_(vc), queue_(container.list_particles())
    {
        shuffle(rng_, queue_);
    }

    bool operator()()
    {
        SGFRD_SCOPE(ns, BDPropagator, this->vc_.access_tracer())
        if(queue_.empty()){return false;}

        // make copy of the next particle
        ParticleID pid; Particle p;
        boost::tie(pid, p) = queue_.back(); queue_.pop_back();
        FaceID fid = this->container_.get_face_id(pid);

        // to restore the position, copy previous state.
        const Real3  prev_pos(p.position());
        const FaceID prev_fid(fid);

        if(this->attempt_reaction(pid, p, fid)){return true;}
        if(p.D() == 0.0)                       {return true;}

        // no 1st kind reaction occured & particle is movable.
        BOOST_AUTO(position,     std::make_pair(p.position(), fid));
        BOOST_AUTO(displacement, draw_displacement(p, fid));
        this->propagate(position, displacement);

        // check escapement and clear volume if needed
        {
            // update local copy of particle
            boost::tie(p.position(), fid) = position;
            if(false == clear_volume(p, fid, pid))
            {
                // rejected. restore position. previous position does not cause
                // overlap because position and species are kept intact.
                ++(this->rejected_move_count_);
                p.position() = prev_pos;
                fid          = prev_fid;
            }
        }

        // retrieve possible reactants (within r1+r2+reaction_length)
        BOOST_AUTO(overlapped, this->list_reaction_overlap(pid, p, fid));

        // check core-overlap
        std::pair<ParticleID, Particle> pp; Real d;
        BOOST_FOREACH(boost::tie(pp, d), overlapped)
        {
            if(d < p.radius() + pp.second.radius())
            {
                // core overlap!
                // restore position and re-collect overlapped particles
                p.position() = prev_pos;
                fid          = prev_fid;
                overlapped   = this->list_reaction_overlap(pid, p, fid);
                break;
            }
        }

        if(overlapped.empty())
        {
            // no reaction-partner exists. overlaps are already cleared. update.
            this->container_.update_particle(pid, p, fid);
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

        BOOST_AUTO(const& rules, this->model_.query_reaction_rules(p.species()));
        SGFRD_TRACE(this->vc_.access_tracer().write(
                    "%1% rules found for particle %2%", rules.size(), pid))
        if(rules.empty()){return false;}

        const Real rnd(this->rng_.uniform(0., 1.));
        SGFRD_TRACE(this->vc_.access_tracer().write(
                    "drawn probability = %1%", rnd))
        Real prob = 0.;
        BOOST_FOREACH(reaction_rule_type const& rule, rules)
        {
            SGFRD_TRACE(this->vc_.access_tracer().write("k * dt = %1%",
                        rule.k() * dt_))
            if((prob += rule.k() * dt_) <= rnd){continue;}
            if(prob >= 1.){std::cerr << "reaction prob exceeds 1" << std::endl;}

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
        BOOST_STATIC_ASSERT(boost::is_same<
            typename boost::iterator_value<Iterator>::type,
            std::pair<std::pair<ParticleID, Particle>, Real> >::value);

        const Real rnd(rng_.uniform(0., 1.));
        Real acc_prob = 0.;

        for(Iterator iter(first); iter != last; ++iter)
        {
            const ParticleID& pid2 = iter->first.first;
            const Particle&   p2   = iter->first.second;

            BOOST_AUTO(const& rules,
                this->model_.query_reaction_rules(p1.species(), p2.species()));
            if(rules.empty())
            {
                // no reaction can occur because there is no rule.
                continue;
            }

            Real acc_prob_local_increase(0.0);
            const Real acceptance_coef = calc_acceptance_coef(p1, p2);
            BOOST_FOREACH(reaction_rule_type const& rule, rules)
            {
                const Real inc = rule.k() * acceptance_coef;
                acc_prob_local_increase += inc;

                acc_prob += inc;
                if(acc_prob <= rnd){continue;}
                if(acc_prob > 1.0)
                {
                    std::cerr << "WARNING: reaction probability exceeds 1\n";
                }

                switch(rule.products().size())
                {
                    case 0: // 2->0 reaction
                    {
                        remove_particle(pid1);
                        remove_particle(pid2);
                        last_reactions_.push_back(std::make_pair(rule,
                             init_reaction_info(pid1, p1, pid2, p2)));
                        return true;
                    }
                    case 1:
                    {
                        const FaceID& f2 =
                            this->container_.get_face_id(pid2);
                        const bool reacted = attempt_reaction_2_to_1(
                            pid1, p1, f1, pid2, p2, f2, std::make_pair(rule,
                                init_reaction_info(pid1, p1, pid2, p2)));
                        if(reacted)
                        {
                            return true;
                        }
                        else
                        {
                            // try next partner
                            // after restoring acceptance probability
                            acc_prob -= acc_prob_local_increase;
                            break;
                        }
                    }
                    default:
                    {
                        throw NotSupported("BDPropagator: 2 -> N (N>1) "
                                           "reaction is not allowed");
                    }
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
        const molecule_info_type mol_info =
            this->container_.get_molecule_info(species_new);
        const Real radius_new = mol_info.radius;
        const Real D_new      = mol_info.D;

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

        const molecule_info_type mol1 = container_.get_molecule_info(sp1);
        const molecule_info_type mol2 = container_.get_molecule_info(sp2);

        const Real D1(mol1.D),      D2(mol2.D),      D12(mol1.D + mol2.D);
        const Real r1(mol1.radius), r2(mol2.radius), r12(mol1.radius + mol2.radius);
        const Real3 n = polygon_.triangle_at(fid).normal();

        if(D1 == 0. && D2 == 0)
        {
            throw NotSupported("BDPropagator::1->2: "
                    "reaction between immobile particles");
        }

        boost::array<std::pair<Real3, FaceID>, 2> newpfs;
        newpfs[0] = std::make_pair(p.position(), fid);
        newpfs[1] = std::make_pair(p.position(), fid);

        {
            const Real3 ipv(draw_ipv(r12, D12, n));
            Real3 disp1(ipv * ( D1 / D12)), disp2(ipv * (-D2 / D12));

            // put two particles next to each other
            this->propagate(newpfs[0], disp1);
            this->propagate(newpfs[1], disp2);

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
        BOOST_AUTO(pp2,
            this->container_.new_particle(particles_new[1], newpfs[1].second));
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
        const species_type     sp_new(rlog.first.products().front());
        const molecule_info_type info(container_.get_molecule_info(sp_new));
        const Real radius_new(info.radius), D_new(info.D);

        const Real3 pos1(p1.position()), pos2(p2.position());
        const Real D1(p1.D()), D2(p2.D());
        const Real D12(D1 + D2);

        // this calculates the center position weighted by Diffusion coef.
        Real3 dp = ecell4::polygon::direction(this->polygon_,
            std::make_pair(pos1, fid1), std::make_pair(pos2, fid2)) * (D1 / D12);

        std::pair<Real3, FaceID> pf1(pos1, fid1);
        this->propagate(pf1, dp);

        if(is_overlapping(pf1, radius_new, pid1, pid2)){return false;}

        const Particle particle_new(sp_new, pf1.first, radius_new, D_new);
        if(false == clear_volume(particle_new, pf1.second, pid1, pid2))
        {
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

             if(std::abs(tilt - M_PI) < 1e-10) return rnd;
        else if(std::abs(tilt + M_PI) < 1e-10) return rnd * (-1.0);

        const Real3 axis = cross_product(Real3(0., 0., 1.), normal);
        return rotate(tilt, axis * (1. / length(axis)), rnd);
    }

    Real3 draw_displacement(const Particle& p, const FaceID& fid)
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
