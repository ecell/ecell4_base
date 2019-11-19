#ifndef ECELL4_SGFRD_SIMULATOR
#define ECELL4_SGFRD_SIMULATOR

#include <greens_functions/GreensFunction2DAbsSym.hpp>
#include <greens_functions/GreensFunction2DRefWedgeAbs.hpp>
#include <greens_functions/GreensFunction2DRadAbs.hpp>

#include <ecell4/core/SimulatorBase.hpp>
#include <ecell4/core/ReactionRule.hpp>
#include <ecell4/core/SerialIDGenerator.hpp>
#include <ecell4/core/geometry.hpp>

#include <ecell4/sgfrd/ShellContainer.hpp>
#include <ecell4/sgfrd/ShellVisitorApplier.hpp>
#include <ecell4/sgfrd/ShellVisitors.hpp>
#include <ecell4/sgfrd/ReactionInfo.hpp>
#include <ecell4/sgfrd/SGFRDEvent.hpp>
#include <ecell4/sgfrd/expected.hpp>

#include <ecell4/sgfrd/tracer.hpp>
#include <ecell4/sgfrd/statistics.hpp>

#include <boost/make_shared.hpp>
#include <boost/container/static_vector.hpp>
#include <boost/container/small_vector.hpp>
#include <tuple>
#include <array>

#include <iostream>

#define SGFRD_LOG(x, y) /**/

namespace ecell4
{
namespace sgfrd
{

class SGFRDSimulator :
    public ecell4::SimulatorBase<SGFRDWorld, ecell4::Model>
{
  public:
    typedef SGFRDSimulator self_type;

    // polygon
    typedef ecell4::Polygon  polygon_type;
    typedef ecell4::Triangle Triangle;
    typedef polygon_type::FaceID   FaceID;
    typedef polygon_type::EdgeID   EdgeID;
    typedef polygon_type::VertexID VertexID;

    // Event & Domain
    typedef SGFRDEvent                      event_type;
    typedef SGFRDEvent::domain_type         domain_type;
    typedef EventID                         event_id_type;
    typedef DomainID                        domain_id_type;
    typedef SGFRDEventScheduler             scheduler_type;
    typedef SGFRDEventScheduler::value_type event_id_pair_type;

    // Simulator
    typedef ecell4::SimulatorBase<SGFRDWorld, ecell4::Model> base_type;
    typedef base_type::world_type world_type;
    typedef base_type::model_type model_type;
    typedef std::pair<ParticleID, Particle> particle_id_pair_type;
    typedef std::tuple<ParticleID, Particle, FaceID> pid_p_fid_tuple_type;
    // length of this type may be tuned
    typedef boost::container::small_vector<pid_p_fid_tuple_type, 5> bursted_type;

    // ShellContainer
    typedef ecell4::SerialIDGenerator<ShellID> shell_id_generator_type;
    typedef ShellContainer                             shell_container_type;
    typedef shell_container_type::shell_type           shell_type;
    typedef shell_container_type::shell_id_pair_type   shell_id_pair_type;
    typedef shell_container_type::circle_type          circle_type;
    typedef shell_container_type::conical_surface_type conical_surface_type;
    typedef shell_container_type::circular_shell_type  circular_shell_type;
    typedef shell_container_type::conical_surface_shell_type
            conical_surface_shell_type;
    typedef shell_visitor_applier<shell_container_type>
            mutable_shell_visitor_applier_type;
    typedef shell_visitor_applier<const shell_container_type>
            immutable_shell_visitor_applier_type;

    // reaction
    typedef ecell4::ReactionRule           reaction_rule_type;
    typedef ReactionInfo                   reaction_info_type;
    typedef std::pair<reaction_rule_type, reaction_info_type> reaction_log_type;
    typedef std::vector<reaction_log_type> reaction_archive_type;

  public:

    SGFRDSimulator(const boost::shared_ptr<world_type>& world,
                   const boost::shared_ptr<model_type>& model,
                   Real bd_dt_factor = 0.01, Real reaction_length = 0.1,
                   const std::string& trace_fname = "sgfrd_trace.log")
        : base_type(world, model), is_dirty_(true), dt_(0),
          bd_dt_factor_(bd_dt_factor), reaction_length_(reaction_length),
          rng_(*(world->rng())), shell_container_(world->polygon()),
          mut_sh_vis_applier(shell_container_), imm_sh_vis_applier(shell_container_),
          tracer_(trace_fname)
    {
        if(1.0 + reaction_length >= single_circular_shell_factor)
        {
            // If this condition is satisfied, sGFRD allows a single shell that
            // has a reactive range outside the shell. But single domain does
            // not consider reactions outside the shell, so it causes an error.
            throw std::invalid_argument("SGFRDSimulator: too large reaction length");
        }
#ifndef ECELL4_SGFRD_NO_TRACE
        tracer_.clear(); // create a file and clear the file if already exists
#endif
    }

    SGFRDSimulator(boost::shared_ptr<world_type> world,
                   Real bd_dt_factor = 0.01, Real reaction_length = 0.1,
                   const std::string& trace_fname = "sgfrd_trace.log")
        : base_type(world), is_dirty_(true), dt_(0),
          bd_dt_factor_(bd_dt_factor), reaction_length_(reaction_length),
          rng_(*(world->rng())), shell_container_(world->polygon()),
          mut_sh_vis_applier(shell_container_), imm_sh_vis_applier(shell_container_),
          tracer_(trace_fname)
    {
        if(1.0 + reaction_length >= single_circular_shell_factor)
        {
            // If this condition is satisfied, sGFRD allows a single shell that
            // has a reactive range outside the shell. But single domain does
            // not consider reactions outside the shell, so it causes an error.
            throw std::invalid_argument("SGFRDSimulator: too large reaction length");
        }
#ifndef ECELL4_SGFRD_NO_TRACE
        tracer_.clear(); // create a file and clear the file if already exists
#endif
    }
    ~SGFRDSimulator() override = default;

    void initialize() override
    {
        if(!(this->is_dirty_))
        {
            // already initialized. do nothing.
            return;
        }
        ParticleID pid; Particle p;
        for(const auto& pidp : this->world_->list_particles())
        {
            std::tie(pid, p) = pidp;
            add_event(create_tight_domain(create_tight_shell(
                      pid, p, this->get_face_id(pid)), pid, p));
        }
        this->is_dirty_ = false;
        return ;
    }
    void finalize()
    {
        assert(this->diagnosis());

        this->is_dirty_ = true;

        const Real tm(this->time());
        while(this->scheduler_.size() != 0)
        {
            this->burst_event(this->scheduler_.pop(), tm);
        }
        return ;
    }
    void finalize(const Real t)
    {
        assert(this->diagnosis());

        this->is_dirty_ = true;

        assert(t < this->next_event_time());
        const Real tm(t);
        this->set_t(t);
        while(this->scheduler_.size() != 0)
        {
            this->burst_event(this->scheduler_.pop(), tm);
        }
        return ;
    }

    Real next_event_time() const
    {
        return this->scheduler_.next_time();
    }

    void step() override
    {
        SGFRD_SCOPE(us, step, tracer_);

        if(this->is_dirty_) // unlikely
        {
            // if this simulator has been finalized (in the last step(upto)),
            // re-initialize it.
            this->initialize();
        }

        this->set_t(this->scheduler_.next_time());
        SGFRD_TRACE(tracer_.write("now t = %1%", this->time()))

        // fire event executes `create_event` inside.
        this->fire_event(this->scheduler_.pop());
        SGFRD_TRACE(tracer_.write("now %1% shells exists", shell_container_.num_shells()))
        SGFRD_TRACE(tracer_.write("now %1% events exists", scheduler_.size()))

        // Set the next dt_.
        // This is required to run this simulator with a FixedIntervalXXXObserver.
        // 1. SimulatorBase::run(obs, upto) checks whether obs.next_time() <
        //    `Simulator::next_time()`.
        // 2. Simulator::next_time() is not virtual. So it cannot be overridden.
        // 3. Simulator::next_time() returns `simulator::t() + simulator::dt()`.
        this->dt_ = this->scheduler_.next_time() - this->time();

        //XXX
        assert(this->diagnosis());

        return;
    }
    /**
     * step and return true if the next time is less than upto.
     * if not, step till upto and return false.
     * @return if the simulator does not rearch upto
     */
    bool step(const Real& upto) override
    {
        if(this->is_dirty_)
        {
            // if this simulator has been finalized (in the last step(upto)),
            // re-initialize it.
            this->initialize();
        }

        if(this->time() > upto)
        {
            // it's too late to stop...
            return false;
        }

        if(this->scheduler_.next_time() < upto)
        {
            // step does re-initialize inside it.
            this->step();
            assert(this->time() < upto);
            return true;
        }
        else if(this->scheduler_.next_time() == upto) // really unlikely
        {
            this->step();
            assert(this->time() == upto);
            this->finalize();
            return false;
        }
        // [[assert axiom: next_time > upto]]

        // burst all the domains at t == upto.
        this->finalize(upto);
        return false;
    }

    void set_t(const Real t) {return this->world_->set_t(t);}

    world_type const& world() const {return *(this->world_);}

    Real dt() const override {return dt_;}
    Real reaction_length() const {return reaction_length_;}

    bool check_reaction() const {return last_reactions_.size() > 0;}
    std::vector<std::pair<ReactionRule, reaction_info_type> > const&
    last_reactions() const {return last_reactions_;}

    bool diagnosis() const;

  private:

    // simple wrappers to call member's member-method (e.g. world_->t()) {{{
    Real uniform_real(){return this->rng_.random();}

    polygon_type const& polygon() const {return *(this->world_->polygon());}

    bool update_particle(const ParticleID& pid, const Particle& p,
                         const FaceID& fid)
    {
        const bool result = this->world_->update_particle(pid, p, fid);
        SGFRD_TRACE(tracer_.write("  particle %1% is updated, %2%", pid, result))
        assert(result == false);
        return result;
    }
    std::pair<std::pair<ParticleID, Particle>, bool>
    create_particle(const Particle& p, const FaceID& fid)
    {
        const std::pair<std::pair<ParticleID, Particle>, bool> result =
            this->world_->new_particle(p, fid);
        assert(result.second);
//         {
//             SGFRD_TRACE(tracer_.write("  particle %1% has been created",
//                         result.first.first))
//         }
//         else
//         {
//             SGFRD_TRACE(tracer_.write("  failed to create particle %1%",
//                         result.first.first))
//         }
        return result;
    }
    void remove_particle(const ParticleID& pid, const FaceID& fid)
    {
        SGFRD_TRACE(tracer_.write("  removing particle %1% on face %2%", pid, fid));
        return this->world_->remove_particle(pid, fid);
    }

    shell_type&       get_shell(ShellID const& id)
    {
        SGFRD_TRACE(tracer_.write("  searching shell %1%", id));
        return shell_container_.get_shell(id);
    }
    shell_type const& get_shell(ShellID const& id) const
    {
        SGFRD_TRACE(tracer_.write("  searching shell %1%", id));
        return shell_container_.get_shell(id);
    }
    void remove_shell(ShellID const& id)
    {
        SGFRD_TRACE(tracer_.write("  removing shell %1%", id));
        return shell_container_.remove_shell(id);
    }
    template<typename shT, typename stridT>
    void update_shell(ShellID const& shid, shT const& sh, stridT strid)
    {
        this->shell_container_.update_shell(shid, sh, strid);
        return;
    }

    FaceID get_face_id(const ParticleID& pid) const
    {
        if(!this->world_->is_on_face(pid))
        {
            SGFRD_TRACE(tracer_.write("  particle %1% is not on face", pid));
            SGFRD_TRACE(tracer_.write("  is particle %1% is in world?: %2%",
                        pid, this->world_->has_particle(pid)));
        }
        return this->world_->get_face_id(pid);
    }
    std::pair<ParticleID, Particle> get_particle(const ParticleID& pid) const
    {return this->world_->get_particle(pid);}

    Real time() const {return this->world_->t();}

    bool event_exists(const event_id_type& id)
    {
        SGFRD_TRACE(tracer_.write("  checking event %1% exists or not", id));
        try
        {
            boost::shared_ptr<event_type> ev = scheduler_.get(id);
            return static_cast<bool>(ev);
        }
        catch(std::out_of_range const& oor)
        {
            return false;
        }
    }

    boost::shared_ptr<event_type> get_event(const event_id_type& id)
    {
        SGFRD_TRACE(tracer_.write("  getting event %1%", id));
        return scheduler_.get(id);
    }
    void remove_event(const event_id_type id)
    {
        SGFRD_TRACE(tracer_.write("  removing event %1%", id));
        this->scheduler_.remove(id);
        return;
    }

    DomainID get_domain_id(Single const& dom) const
    {
        return boost::apply_visitor(domain_id_getter(), get_shell(dom.shell_id()));
    }
    DomainID get_domain_id(Pair const& dom) const
    {
        return boost::apply_visitor(domain_id_getter(), get_shell(dom.shell_id()));
    }
    DomainID get_domain_id(Multi const& dom) const
    {
        return boost::apply_visitor(domain_id_getter(),
                                    get_shell(dom.shell_ids().front()));
    }
    // }}}

  private:

//----------------------------------- single -----------------------------------

    /*! execute Event associated with Domain, remove shell, create next event. */
    void          fire_single(const Single& dom, DomainID did);
    bursted_type burst_single(const Single& dom, const Real tm);

    std::tuple<ParticleID, Particle, FaceID>
    propagate_single(
            const shell_type& sh, const Single& dom, const Real tm)
    {
        switch(sh.which())
        {
            case shell_container_type::circular_shell:
            {
                return propagate_single_circular(
                    boost::get<circular_shell_type>(sh), dom, tm);
            }
            case shell_container_type::conical_shell:
            {
                return propagate_single_conical(
                    boost::get<conical_surface_shell_type>(sh), dom, tm);
            }
            default:
            {
                throw std::logic_error((boost::format(
                    "boost::variant<shells>::which(): invalid value(%1%)") %
                    sh.which()).str());
            }
        }
    }

    std::tuple<ParticleID, Particle, FaceID>
    propagate_single_circular(
        const circular_shell_type& sh, const Single& dom, const Real tm);
    std::tuple<ParticleID, Particle, FaceID>
    propagate_single_conical(
        const conical_surface_shell_type& sh, const Single& dom, const Real tm);

    std::tuple<ParticleID, Particle, FaceID>
    escape_single(const shell_type& sh, const Single& dom)
    {
        switch(sh.which())
        {
            case shell_container_type::circular_shell:
            {
                return escape_single_circular(
                    boost::get<circular_shell_type>(sh), dom);
            }
            case shell_container_type::conical_shell:
            {
                return escape_single_conical(
                    boost::get<conical_surface_shell_type>(sh), dom);
            }
            default:
            {
                throw std::logic_error((boost::format(
                    "boost::variant<shells>::which(): invalid value(%1%)") %
                    sh.which()).str());
            }
        }
    }

    std::tuple<ParticleID, Particle, FaceID>
    escape_single_circular(const circular_shell_type& sh, const Single& dom);
    std::tuple<ParticleID, Particle, FaceID>
    escape_single_conical(const conical_surface_shell_type& sh, const Single& dom);

    boost::container::static_vector<pid_p_fid_tuple_type, 2>
    reaction_single(const shell_type& sh, const Single& dom, const DomainID did);

    boost::container::static_vector<pid_p_fid_tuple_type, 2>
    attempt_reaction_single(
            const shell_type& sh,  const DomainID did, const Single& dom,
            const ParticleID& pid, const Particle&  p, const FaceID& fid);

    boost::container::static_vector<pid_p_fid_tuple_type, 2>
    attempt_reaction_1_to_1(const ReactionRule& rule,
            const shell_type& sh,  const DomainID did, const Single& dom,
            const ParticleID& pid, const Particle&  p, const FaceID& fid);

    boost::container::static_vector<pid_p_fid_tuple_type, 2>
    attempt_reaction_1_to_2(const ReactionRule& rule,
            const shell_type& sh,  const DomainID did, const Single& dom,
            const ParticleID& pid, const Particle&  p, const FaceID& fid);

    std::pair<ShellID, circle_type>
    create_single_circular_shell(
            const std::pair<Real3, FaceID>& pos, const Real size,
            bool check = true)
    {
        SGFRD_SCOPE(ns, create_single_circular_shell, tracer_);
        SGFRD_TRACE(tracer_.write("shell size = %1%", size))

        const ShellID id(shell_id_gen());
        const circle_type shape(size, pos.first,
                                polygon().triangle_at(pos.second).normal());
        if(check)
        {
            shell_container_.check_add_shell(
                id, circular_shell_type(shape, pos.second), pos.second,
                "create_single_circular_shell");
        }
        else
        {
             shell_container_.add_shell(
                id, circular_shell_type(shape, pos.second), pos.second);
        }
        SGFRD_TRACE(tracer_.write("the shell id is %1%", id))
        return std::make_pair(id, shape);
    }
    std::pair<ShellID, conical_surface_type>
    create_single_conical_surface_shell(
            const VertexID& vid, const Real size, bool check = true)
    {
        SGFRD_SCOPE(ns, create_single_conical_shell, tracer_);
        SGFRD_TRACE(tracer_.write("vertex ID  = %1%", vid));
        SGFRD_TRACE(tracer_.write("shell size = %1%", size));

        const ShellID id(shell_id_gen());
        const conical_surface_type shape(polygon().position_at(vid),
                                         polygon().apex_angle_at(vid), size);
        if(check)
        {
            shell_container_.check_add_shell(
                    id, conical_surface_shell_type(shape, vid), vid,
                    "create_single_conical_shell");
        }
        else // no check required.
        {
            shell_container_.add_shell(
                    id, conical_surface_shell_type(shape, vid), vid);
        }
        SGFRD_TRACE(tracer_.write("the shell id is %1%", id));

        return std::make_pair(id, shape);
    }

    Single create_single(const std::pair<ShellID, circle_type>& sh,
                         const ParticleID& pid, const Particle& p)
    {
        SGFRD_SCOPE(ns, create_single_circular_domain, tracer_);
        SGFRD_TRACE(tracer_.write("shell shape = %1%", sh.second))

        const greens_functions::GreensFunction2DAbsSym
            gf(/* D = */ p.D(),
               /* a = */ sh.second.size() - p.radius());
        const Real t_escape = gf.drawTime(uniform_real());
        SGFRD_TRACE(tracer_.write("calculated escape_time = %1%", t_escape))

        //consider whether reaction occur or not.
        const Real t_reaction = draw_reaction_time(calc_k_tot(
            this->model_->query_reaction_rules(p.species())));
        SGFRD_TRACE(tracer_.write("calculated reaction_time = %1%", t_reaction))

        if(t_reaction < t_escape)
        {
            SGFRD_TRACE(tracer_.write("single event is set as reaction"))
            return Single(Single::REACTION, t_reaction, this->time(), sh.first,
                          std::make_pair(pid, p));
        }
        else
        {
            SGFRD_TRACE(tracer_.write("single event is set as escape"))
            return Single(Single::ESCAPE, t_escape, this->time(), sh.first,
                          std::make_pair(pid, p));
        }
    }

    Single create_single(const std::pair<ShellID, conical_surface_type>& sh,
                         const ParticleID& pid, const Particle& p)
    {
        SGFRD_SCOPE(ns, create_single_conical_domain, tracer_);
        SGFRD_TRACE(tracer_.write("shell shape = %1%", sh.second))
        SGFRD_TRACE(tracer_.write("arguments: shell ID = %2%, shell size = %2%"
            ", particle ID = %3%, Particle position = %4%",
            sh.first, sh.second.size(), pid, p.position()))

        const Real D   = p.D();
        const Real a   = sh.second.size() - p.radius();
        const Real phi = sh.second.apex_angle();
        const Real r0  = length(
                this->polygon().periodic_transpose(p.position(), sh.second.apex()) -
                sh.second.apex());

        SGFRD_TRACE(tracer_.write("D   = %1%", D  ))
        SGFRD_TRACE(tracer_.write("r0  = %1%", r0 ))
        SGFRD_TRACE(tracer_.write("a   = %1%", a  ))
        SGFRD_TRACE(tracer_.write("phi = %1%", phi))

        const greens_functions::GreensFunction2DRefWedgeAbs gf(D, r0, a, phi);
        const Real t_escape = gf.drawTime(uniform_real());

        SGFRD_TRACE(tracer_.write("calculated escape_time = %1%", t_escape))

        const Real t_reaction = draw_reaction_time(calc_k_tot(
            this->model_->query_reaction_rules(p.species())));
        SGFRD_TRACE(tracer_.write("calculated reaction_time = %1%", t_reaction))

        if(t_reaction < t_escape)
        {
            SGFRD_TRACE(tracer_.write("single event is set as reaction"))
            return Single(Single::REACTION, t_reaction, this->time(), sh.first,
                          std::make_pair(pid, p));
        }
        else
        {
            SGFRD_TRACE(tracer_.write("single event is set as escape"))
            return Single(Single::ESCAPE, t_escape, this->time(), sh.first,
                          std::make_pair(pid, p));
        }
    }

//------------------------------------ pair ------------------------------------

    void fire_pair(Pair& dom, DomainID did)
    {
        SGFRD_SCOPE(us, fire_pair, tracer_);
        const greens_functions::GreensFunction2DRadAbs
            gf_ipv(dom.D_ipv(), dom.kf(), dom.r0(), dom.sigma(), dom.R_ipv());
        if(dom.eventkind() == Pair::IV_UNDETERMINED)
        {
            SGFRD_TRACE(tracer_.write("pair event is IV_UNDERTERMINED. "
                        "determine iv event kind here."))
            const greens_functions::GreensFunction::EventKind iv_kind =
                gf_ipv.drawEventType(this->uniform_real(), dom.dt());
            if(iv_kind == greens_functions::GreensFunction::IV_ESCAPE)
            {
                dom.eventkind() = Pair::IV_ESCAPE;
                SGFRD_TRACE(tracer_.write("pair event kind = IV_ESCAPE"));
            }
            else if(iv_kind == greens_functions::GreensFunction::IV_REACTION)
            {
                dom.eventkind() = Pair::IV_REACTION;
                SGFRD_TRACE(tracer_.write("pair event kind = IV_REACTION"));
            }
            else
            {
                throw std::runtime_error("neither Pair::IV_ESCAPE/REACTION");
            }
        }

        boost::optional<std::size_t> reactant_index(boost::none);
        const ShellID sid(dom.shell_id());
        switch(dom.eventkind())
        {
            case Pair::SINGLE_REACTION_1:
            {
                SGFRD_SCOPE(ns, case_SINGLE_REACTION_1, tracer_);
                SGFRD_TRACE(tracer_.write("pair event kind = SINGLE_REACTION_1"));
                reactant_index = 0;
                // DO NOT BREAK SWITCH-CASE HERE!
            }
            case Pair::SINGLE_REACTION_2:
            {
                SGFRD_SCOPE(ns, case_SINGLE_REACTION_2, tracer_);
                if(!reactant_index)
                {
                    SGFRD_TRACE(tracer_.write(
                                "pair event kind = SINGLE_REACTION_2"));
                    reactant_index = 1;
                }

                // first, update 2 particles
                std::array<std::tuple<ParticleID, Particle, FaceID>, 2>
                    propagated = this->propagate_pair(
                        this->get_shell(sid), dom, this->time());
                SGFRD_TRACE(tracer_.write("particles are propagated"));

                this->remove_shell(sid);
                SGFRD_TRACE(tracer_.write("shell %1% removed", sid));

                std::array<ShellID,    2> sids;
                std::array<DomainID,   2> dids;
                std::array<Single,     2> doms;
                std::array<ParticleID, 2> pids;
                std::array<Particle,   2> ps;
                std::array<FaceID,     2> fids;

                // add tight-domain for them to detect overlap
                for(std::size_t i=0; i<2; ++i)
                {
                    std::tie(pids[i], ps[i], fids[i]) = propagated[i];
                    SGFRD_TRACE(tracer_.write("adding tight domain > %1%", pids[i]))
                    sids[i] = create_tight_shell(pids[i], ps[i], fids[i]);
                    doms[i] = create_tight_domain(sids[i], pids[i], ps[i]);
                    dids[i] = add_event(doms[i]);
                }
                SGFRD_TRACE(tracer_.write("tight-domains assigned"));

                // after that, attempt single reaction
                const std::size_t ridx = *reactant_index;
                SGFRD_TRACE(tracer_.write("reactant_index = %1%", ridx));

                auto results = this->attempt_reaction_single(
                           this->get_shell(sids[ridx]), dids[ridx], doms[ridx],
                           pids[ridx], ps[ridx], fids[ridx]);

                if(results.size() != 1 || std::get<0>(results.front()) != pids[ridx])
                {
                    STAT(stat_reaction_condition.add_count(PairFirstOrder));
                }
                else // reaction fails.
                {
                    STAT(stat_reaction_condition.add_count(PairFirstOrderFailed));
                }

                this->remove_shell(sids[ridx]);
                SGFRD_TRACE(tracer_.write("shell %1% removed", sids[ridx]));
                this->remove_event(dids[ridx]);
                SGFRD_TRACE(tracer_.write("event %1% removed", dids[ridx]));

                SGFRD_TRACE(tracer_.write("reaction attempted"));

                // add domain to each reactant
                ParticleID pid; Particle p; FaceID fid;
                for(const auto& pidpf : results)
                {
                    std::tie(pid, p, fid) = pidpf;
                    SGFRD_TRACE(tracer_.write("adding next event for %1%", pid))
                    add_event(create_tight_domain(
                              create_tight_shell(pid, p, fid), pid, p));
                }
                return;
            }
            case Pair::COM_ESCAPE: // fallthrough
            case Pair::IV_ESCAPE:
            {
                SGFRD_SCOPE(ns, case_COM_OR_IV_ESCAPE, tracer_);
                std::array<std::tuple<ParticleID, Particle, FaceID>, 2>
                    escaped = this->escape_pair( // dispatch com/ipv here
                        this->get_shell(sid), dom, this->time());
                SGFRD_TRACE(tracer_.write("particles escaped"));

                this->remove_shell(sid);
                SGFRD_TRACE(tracer_.write("shell %1% removed", sid));

                ParticleID pid; Particle p; FaceID fid;
                for(const auto& pidpf : escaped)
                {
                    std::tie(pid, p, fid) = pidpf;
                    SGFRD_TRACE(tracer_.write("adding next event for %1%", pid))
                    add_event(create_tight_domain(
                                create_tight_shell(pid, p, fid), pid, p));
                }
                return;
            }
            case Pair::IV_REACTION:
            {
                SGFRD_SCOPE(ns, case_IV_REACTION, tracer_);
                boost::container::small_vector<
                    std::tuple<ParticleID, Particle, FaceID>, 2>
                        products = this->attempt_pair_reaction(
                                this->get_shell(sid), dom, this->time());
                SGFRD_TRACE(tracer_.write("iv reaction attempted"));
                this->remove_shell(sid);
                SGFRD_TRACE(tracer_.write("shell %1% removed", sid));

                if(!products.empty())
                {
                    ParticleID pid; Particle p; FaceID fid;
                    for(const auto& pidpf : products)
                    {
                        std::tie(pid, p, fid) = pidpf;
                        SGFRD_TRACE(tracer_.write("adding next event for %1%", pid));
                        add_event(create_tight_domain(
                                create_tight_shell(pid, p, fid), pid, p));
                    }
                }
                return;
            }
            default:
            {
                throw std::invalid_argument("fire_pair(): invalid event kind");
            }
        }
    }

    bursted_type burst_pair(const Pair& dom, const Real tm)
    {
        SGFRD_SCOPE(us, burst_pair, tracer_);
        const ShellID sid(dom.shell_id());
        SGFRD_TRACE(tracer_.write("pair shell id = %1%", sid));

        // no reaction occurs because `burst` occurs before any other event.
        bursted_type results;
        std::array<std::tuple<ParticleID, Particle, FaceID>, 2>
            propagated(this->propagate_pair(this->get_shell(sid), dom, tm));
        results.push_back(propagated[0]);
        results.push_back(propagated[1]);

        SGFRD_TRACE(tracer_.write("particle %1% and %2% propagated",
            std::get<0>(propagated[0]), std::get<0>(propagated[1])));

        this->remove_shell(sid);
        SGFRD_TRACE(tracer_.write("shell(%1%) removed", sid));
        return results;
    }

    std::array<std::tuple<ParticleID, Particle, FaceID>, 2>
    propagate_pair(const shell_type& sh, const Pair& dom, const Real tm)
    {
        SGFRD_SCOPE(us, propagate_pair, tracer_);
        switch(sh.which())
        {
            case shell_container_type::circular_shell:
            {
                return propagate_circular_pair(
                        boost::get<circular_shell_type>(sh), dom, tm);
            }
            default:
            {
                throw std::logic_error((boost::format(
                    "propagate_pair: shell::which() returns invalid value(%1%)")
                    % sh.which()).str());
            }
        }
    }

    std::array<std::tuple<ParticleID, Particle, FaceID>, 2>
    propagate_circular_pair(const circular_shell_type& sh,
            const Pair& dom, const Real tm)
    {
        SGFRD_SCOPE(us, propagate_circular_pair, tracer_);

        // calculate displacements
        const Real dt = tm - dom.begin_time();
        const greens_functions::GreensFunction2DRadAbs
            gf_ipv(dom.D_ipv(), dom.kf(), dom.r0(), dom.sigma(), dom.R_ipv());
        const Real l_ipv     = gf_ipv.drawR(this->uniform_real(), dt);
        const Real theta_ipv = (this->uniform_real() < 0.5 ? -1.0 : 1.0) *
                               gf_ipv.drawTheta(this->uniform_real(), l_ipv, dt);

        const greens_functions::GreensFunction2DAbsSym
            gf_com(dom.D_com(), dom.R_com());
        const Real l_com     = gf_com.drawR(this->uniform_real(), dt);
        const Real theta_com = this->uniform_real() *
                               boost::math::constants::two_pi<Real>();

        const FaceID       sh_fid = sh.structure_id();
        const Triangle&         f = this->polygon().triangle_at(sh_fid);
        const Real3 direction_com = rotate(theta_com, f.normal(), f.represent());

        const Real3 disp_com = direction_com * (l_com / length(direction_com));
        const Real3 disp_ipv = rotate(theta_ipv, f.normal(), dom.ipv()) *
                               (l_ipv / length(dom.ipv()));

        // update position
        // XXX to treat polygon surface structure, calculate displacement for
        // each partile without polygon information first
        // and then apply it to each particle with structure

        Particle         p1   = dom.particle_at(0);
        Particle         p2   = dom.particle_at(1);
        const ParticleID pid1 = dom.particle_id_at(0);
        const ParticleID pid2 = dom.particle_id_at(1);

        // ipv is a vector from p1 to p2
        const Real  rD12        = 1.0 / (p1.D() + p2.D());
        const Real3 disp_ipv_p1 = disp_ipv * (-p1.D() * rD12);
        const Real3 disp_ipv_p2 = disp_ipv * ( p2.D() * rD12);

        Real3 disp_p1 = disp_com + disp_ipv_p1;
        Real3 disp_p2 = disp_com + disp_ipv_p2;
        // start from CoM
        std::pair<Real3, FaceID> pos_p1(sh.position(), sh.structure_id());
        std::pair<Real3, FaceID> pos_p2(sh.position(), sh.structure_id());

        pos_p1 = ecell4::polygon::travel(this->polygon(), pos_p1, disp_p1, 2);
        pos_p2 = ecell4::polygon::travel(this->polygon(), pos_p2, disp_p2, 2);

        p1.position() = pos_p1.first;
        p2.position() = pos_p2.first;

        this->update_particle(pid1, p1, pos_p1.second);
        this->update_particle(pid2, p2, pos_p2.second);

        std::array<std::tuple<ParticleID, Particle, FaceID>, 2> results;
        results[0] = std::make_tuple(pid1, p1, pos_p1.second);
        results[1] = std::make_tuple(pid2, p2, pos_p2.second);

        // check particles are still inside of the shell after this propagation
        assert(
            ecell4::polygon::distance(this->polygon(), pos_p1,
                std::make_pair(sh.position(), sh.structure_id())) <=
            sh.size() - p1.radius()
        );
        assert(
            ecell4::polygon::distance(this->polygon(), pos_p2,
                std::make_pair(sh.position(), sh.structure_id())) <=
            sh.size() - p2.radius()
        );
        return results;
    }

    std::array<std::tuple<ParticleID, Particle, FaceID>, 2>
    escape_pair(const shell_type& sh, const Pair& dom, const Real tm)
    {
        SGFRD_SCOPE(us, escape_pair, tracer_);
        switch(dom.eventkind())
        {
            case Pair::COM_ESCAPE:
            {
                return escape_com_pair(sh, dom, tm);
            }
            case Pair::IV_ESCAPE:
            {
                return escape_ipv_pair(sh, dom, tm);
            }
            default:
            {
                throw std::invalid_argument("escape_pair(): invalid event kind");
            }
        }
    }
    // center of mass escapes from inner-domain
    std::array<std::tuple<ParticleID, Particle, FaceID>, 2>
    escape_com_pair(const shell_type& sh, const Pair& dom, const Real tm)
    {
        SGFRD_SCOPE(us, escape_com_pair, tracer_);
        switch(sh.which())
        {
            case shell_container_type::circular_shell:
            {
                return escape_com_circular_pair(
                    boost::get<circular_shell_type>(sh), dom, tm);
            }
            default:
            {
                throw std::logic_error((boost::format(
                    "boost::variant<shells>::which(): invalid value(%1%)") %
                    sh.which()).str());
            }
        }
    }
    // inter particle vector become too long to react each other
    std::array<std::tuple<ParticleID, Particle, FaceID>, 2>
    escape_com_circular_pair(
            const circular_shell_type& sh, const Pair& dom, const Real tm)
    {
        SGFRD_SCOPE(us, escape_com_circular_pair, tracer_);

        // calculate displacements
        const Real dt = tm - dom.begin_time();
        const greens_functions::GreensFunction2DRadAbs
            gf_ipv(dom.D_ipv(), dom.kf(), dom.r0(), dom.sigma(), dom.R_ipv());

        const Real l_ipv     = gf_ipv.drawR    (this->uniform_real(), dt);
        const Real theta_ipv = gf_ipv.drawTheta(this->uniform_real(), l_ipv, dt) *
                               (uniform_real() < 0.5 ? -1 : 1) ;
        // drawTheta returns value in [0, pi]. we need to determine the directionality.
        const Real l_com     = dom.R_com();
        const Real theta_com = this->uniform_real() *
                               boost::math::constants::two_pi<Real>();

        SGFRD_TRACE(tracer_.write("r_ipv = %1%, r_com = %2%", dom.R_ipv(), dom.R_com()));
        SGFRD_TRACE(tracer_.write("l_ipv = %1%, theta_ipv = %2%", l_ipv, theta_ipv));
        SGFRD_TRACE(tracer_.write("l_com = %1%, theta_com = %2%", l_com, theta_com));

        const FaceID sh_fid = sh.structure_id();
        const Triangle&   f = this->polygon().triangle_at(sh_fid);

        const Real3 direction_com = rotate(theta_com, f.normal(), f.represent());
        const Real3 disp_com = direction_com * (l_com / length(direction_com));
        const Real3 disp_ipv = rotate(theta_ipv, f.normal(), dom.ipv()) *
                               (l_ipv / length(dom.ipv()));
        SGFRD_TRACE(tracer_.write("length of disp_ipv = %1%", length(disp_ipv)));

        // update position
        // XXX to treat polygon surface structure, calculate displacement for
        // each partile without polygon information first
        // and then apply it to each particle with structure

        Particle         p1   = dom.particle_at(0);
        Particle         p2   = dom.particle_at(1);
        const ParticleID pid1 = dom.particle_id_at(0);
        const ParticleID pid2 = dom.particle_id_at(1);

        const std::pair<Real3, FaceID> pos_com(sh.position(), sh.structure_id());

        // ipv is a vector from p1 to p2
        const Real3 disp_p1 = disp_com + disp_ipv * (-p1.D() / (p1.D() + p2.D()));
        const Real3 disp_p2 = disp_com + disp_ipv * ( p2.D() / (p1.D() + p2.D()));
        std::pair<Real3, FaceID> pos_p1(pos_com);
        std::pair<Real3, FaceID> pos_p2(pos_com);
        SGFRD_TRACE(tracer_.write("p1 is on face %1%, p2 is on face %2%",
                                  get_face_id(pid1), get_face_id(pid2)))

        pos_p1 = ecell4::polygon::travel(this->polygon(), pos_p1, disp_p1, 2);
        pos_p2 = ecell4::polygon::travel(this->polygon(), pos_p2, disp_p2, 2);

        // check distance between p1 and p2; it should be the same as l_ipv
        SGFRD_TRACE(tracer_.write("distance after travel = %1%",
                    ecell4::polygon::distance(this->polygon(), pos_p1, pos_p2)))
        assert(std::abs(ecell4::polygon::distance(
                this->polygon(), pos_p1, pos_p2,  true) - l_ipv) < l_ipv * 1e-6);

        p1.position() = pos_p1.first;
        p2.position() = pos_p2.first;

        this->update_particle(pid1, p1, pos_p1.second);
        this->update_particle(pid2, p2, pos_p2.second);

        std::array<std::tuple<ParticleID, Particle, FaceID>, 2> results;
        results[0] = std::make_tuple(pid1, p1, pos_p1.second);
        results[1] = std::make_tuple(pid2, p2, pos_p2.second);

        // check particles are still inside of the shell after this propagation
        assert(
            ecell4::polygon::distance(this->polygon(), pos_p1,
                std::make_pair(sh.position(), sh.structure_id())) <=
            sh.size() - p1.radius()
        );
        assert(
            ecell4::polygon::distance(this->polygon(), pos_p2,
                std::make_pair(sh.position(), sh.structure_id())) <=
            sh.size() - p2.radius()
        );

        return results;
    }

    std::array<std::tuple<ParticleID, Particle, FaceID>, 2>
    escape_ipv_pair(const shell_type& sh, const Pair& dom, const Real tm)
    {
        SGFRD_SCOPE(us, escape_ipv_pair, tracer_);
        switch(sh.which())
        {
            case shell_container_type::circular_shell:
            {
                return escape_ipv_circular_pair(
                    boost::get<circular_shell_type>(sh), dom, tm);
            }
            default:
            {
                throw std::logic_error((boost::format(
                    "boost::variant<shells>::which(): invalid value(%1%)") %
                    sh.which()).str());
            }
        }
    }
    std::array<std::tuple<ParticleID, Particle, FaceID>, 2>
    escape_ipv_circular_pair(
            const circular_shell_type& sh, const Pair& dom, const Real tm)
    {
        SGFRD_SCOPE(us, escape_ipv_circular_pair, tracer_);
        // calculate length and direction of inter particle vector
        const Real dt = tm - dom.begin_time();
        const greens_functions::GreensFunction2DRadAbs
            gf_ipv(dom.D_ipv(), dom.kf(), dom.r0(), dom.sigma(), dom.R_ipv());
        const Real l_ipv     = dom.R_ipv();
        const Real theta_ipv = gf_ipv.drawTheta(this->uniform_real(), l_ipv, dt) *
                               (uniform_real() < 0.5 ? -1 : 1) ;

        // calculate position of the center of mass
        const greens_functions::GreensFunction2DAbsSym
            gf_com(dom.D_com(), dom.R_com());
        const Real l_com     = gf_com.drawR(this->uniform_real(), dt);
        const Real theta_com = this->uniform_real() *
                               boost::math::constants::two_pi<Real>();

        assert(l_com <= dom.R_com());

        SGFRD_TRACE(tracer_.write("r_ipv = %1%, r_com = %2%", dom.R_ipv(), dom.R_com()));
        SGFRD_TRACE(tracer_.write("l_ipv = %1%, theta_ipv = %2%", l_ipv, theta_ipv));
        SGFRD_TRACE(tracer_.write("l_com = %1%, theta_com = %2%", l_com, theta_com));

        const FaceID sh_fid = sh.structure_id();
        const Triangle&   f = this->polygon().triangle_at(sh_fid);
        const Real3 direction_com = rotate(theta_com, f.normal(), f.represent());

        const Real3 disp_com = direction_com * (l_com / length(direction_com));
        const Real3 disp_ipv = rotate(theta_ipv, f.normal(), dom.ipv()) *
                               (l_ipv / length(dom.ipv()));
        SGFRD_TRACE(tracer_.write("length of disp_com = %1%", length(disp_com)));
        SGFRD_TRACE(tracer_.write("length of disp_ipv = %1%", length(disp_ipv)));

        // update position
        // XXX to treat polygon surface structure, calculate displacement for
        // each partile without polygon information first
        // and then apply it to each particle with structure

        Particle         p1   = dom.particle_at(0);
        Particle         p2   = dom.particle_at(1);
        const ParticleID pid1 = dom.particle_id_at(0);
        const ParticleID pid2 = dom.particle_id_at(1);

        assert(!(p1.D() == 0.0 && p2.D() == 0.0));
        SGFRD_TRACE(tracer_.write("previous position of p1 = %1%", p1.position()));
        SGFRD_TRACE(tracer_.write("previous position of p2 = %1%", p2.position()));

        const Real3&  pos_com(sh.position());
        const FaceID& fid_com(sh.structure_id());
        // ipv is a vector from p1 to p2
        Real3 disp_p1 = disp_com + disp_ipv * (-p1.D() / (p1.D() + p2.D()));
        Real3 disp_p2 = disp_com + disp_ipv * ( p2.D() / (p1.D() + p2.D()));
        std::pair<Real3, FaceID> pos_p1(pos_com, fid_com);
        std::pair<Real3, FaceID> pos_p2(pos_com, fid_com);

        pos_p1 = ecell4::polygon::travel(this->polygon(), pos_p1, disp_p1, 2);
        pos_p2 = ecell4::polygon::travel(this->polygon(), pos_p2, disp_p2, 2);

        SGFRD_TRACE(tracer_.write("position of p1 after reaction = %1%", pos_p1.first));
        SGFRD_TRACE(tracer_.write("position of p2 after reaction = %1%", pos_p2.first));

        p1.position() = pos_p1.first;
        p2.position() = pos_p2.first;

        this->update_particle(pid1, p1, pos_p1.second);
        this->update_particle(pid2, p2, pos_p2.second);

        std::array<std::tuple<ParticleID, Particle, FaceID>, 2> results;
        results[0] = std::make_tuple(pid1, p1, pos_p1.second);
        results[1] = std::make_tuple(pid2, p2, pos_p2.second);

        // check particles are still inside of the shell after this propagation
        assert(
            ecell4::polygon::distance(this->polygon(), pos_p1,
                std::make_pair(sh.position(), sh.structure_id())) <=
            sh.size() - p1.radius()
        );
        assert(
            ecell4::polygon::distance(this->polygon(), pos_p2,
                std::make_pair(sh.position(), sh.structure_id())) <=
            sh.size() - p2.radius()
        );

        return results;
    }

    boost::container::small_vector<std::tuple<ParticleID, Particle, FaceID>, 2>
    attempt_pair_reaction(const shell_type& sh, const Pair& dom, const Real tm)
    {
        SGFRD_SCOPE(us, attempt_pair_reaction, tracer_);
        switch(sh.which())
        {
            case shell_container_type::circular_shell:
            {
                return attempt_circular_pair_reaction(
                    boost::get<circular_shell_type>(sh), dom, tm);
            }
            default:
            {
                throw std::logic_error((boost::format(
                    "boost::variant<shells>::which(): invalid value(%1%)") %
                    sh.which()).str());
            }
        }
    }

    boost::container::small_vector<std::tuple<ParticleID, Particle, FaceID>, 2>
    attempt_circular_pair_reaction(
            const circular_shell_type& sh, const Pair& dom, const Real tm)
    {
        SGFRD_SCOPE(us, attempt_circular_pair_reaction, tracer_);
        // 1. update one with com diffusion
        // 2. mutate the particle to product
        // 3. and remove the other one

        Particle p1 = dom.particle_at(0);
        Particle p2 = dom.particle_at(1);
        const ParticleID pid1 = dom.particle_id_at(0);
        const ParticleID pid2 = dom.particle_id_at(1);
        const FaceID fid1 = this->get_face_id(pid1);
        const FaceID fid2 = this->get_face_id(pid2);

        const auto& rules = this->model_->query_reaction_rules(
                    p1.species(), p2.species());
        assert(!rules.empty());

        const Real k_tot = this->calc_k_tot(rules);
        boost::optional<ReactionRule const&> optr(boost::none);
        if(rules.size() != 1)
        {
            Real rndr = this->uniform_real() * k_tot;
            for(const auto& rl : rules)
            {
                rndr -= rl.k();
                if(rndr < 0.0)
                {
                    optr = rl;
                    break;
                }
            }
            // optr maybe empty because of numerical error. in that case,
            // use domain.reactants().back().
            // if domain.reactions().size() == 1, it is okay to use
            // domain.reactants().back() because it is the only rule
            // that can be applied.
        }
        ReactionRule const& rule = static_cast<bool>(optr) ? *optr : (rules.back());

        switch(rule.products().size())
        {
            case 0:
            {
                SGFRD_TRACE(tracer_.write("degradation reaction occurs."))

                this->remove_particle(pid1, fid1);
                this->remove_particle(pid2, fid2);
                last_reactions_.push_back(std::make_pair(rule,
                    make_degradation_reaction_info(tm, pid1, p1, pid2, p2)));

                STAT(stat_reaction_condition.add_count(PairSecondOrder));

                return boost::container::small_vector<
                    std::tuple<ParticleID, Particle, FaceID>, 2>(0ul);
            }
            case 1:
            {
                SGFRD_TRACE(tracer_.write("2->1 reaction occurs."))

                // calculate next position of particle pair
                const Real dt = tm - dom.begin_time();
                const greens_functions::GreensFunction2DAbsSym
                    gf_com(dom.D_com(), dom.R_com());
                const Real l_com     = gf_com.drawR(this->uniform_real(), dt);
                const Real theta_com = this->uniform_real() *
                                       boost::math::constants::two_pi<Real>();

                const FaceID sh_fid = sh.structure_id();
                const Triangle&   f = this->polygon().triangle_at(sh_fid);
                const Real3 direction_com =
                    rotate(theta_com, f.normal(), f.represent());

                Real3 disp_com = direction_com * (l_com / length(direction_com));
                std::pair<Real3, FaceID> pos_com(sh.position(), sh_fid);

                pos_com = ecell4::polygon::travel(this->polygon(), pos_com, disp_com, 2);

                const Real3  pos_new = pos_com.first;
                const FaceID fid_new = pos_com.second;

                // make new particle
                const auto species_new = rule.products().front();
                const auto molinfo     = this->world_->get_molecule_info(species_new);
                const Real radius_new  = molinfo.radius;
                const Real D_new       = molinfo.D;
                const Particle p_new(species_new, pos_new, radius_new, D_new);

                inside_checker is_inside_of(
                    p_new.position(), p_new.radius(), fid_new, this->polygon());
                if(!is_inside_of(sh))
                {
                    // particle sticks out from the shell after reaction
                    // because of its radius.

                    SGFRD_SCOPE(us, particle_goes_outside, tracer_)
                    const DomainID did = this->get_domain_id(dom);

                    const bool no_overlap =
                        this->burst_and_shrink_overlaps(p_new, fid_new, did);
                    SGFRD_TRACE(tracer_.write("no_overlap = %1%", no_overlap))

                    // particle overlaps with other particle/domain outside the
                    // shell. rollback the state.
                    if(!no_overlap)
                    {
                        SGFRD_TRACE(tracer_.write(
                                    "reject the reaction because of no space"))
                        const greens_functions::GreensFunction2DRadAbs
                            gf_ipv(dom.D_ipv(), dom.kf(), dom.r0(), dom.sigma(),
                                   dom.R_ipv());
                        const Real theta_ipv = gf_ipv.drawTheta(
                                this->uniform_real(), dom.sigma(), dt) *
                               (uniform_real() < 0.5 ? -1 : 1);
                        const Real3 disp_ipv =
                            rotate(theta_ipv, f.normal(), dom.ipv());

                        // TODO make function like
                        // propagate_pair(ipv0, disp_com, ipv1)
                        const Real  ratio_p1     = -p1.D() / (p1.D() + p2.D());
                        const Real  ratio_p2     =  p2.D() / (p1.D() + p2.D());
                        const Real3 disp_ipv_p1  = disp_ipv  * ratio_p1;
                        const Real3 disp_ipv_p2  = disp_ipv  * ratio_p2;
                        const Real3 disp_ipv0_p1 = dom.ipv() * ratio_p1;
                        const Real3 disp_ipv0_p2 = dom.ipv() * ratio_p2;

                        Real3 disp_p1 = disp_com - disp_ipv0_p1 + disp_ipv_p1;
                        Real3 disp_p2 = disp_com - disp_ipv0_p2 + disp_ipv_p2;
                        disp_p1 *= (1. + minimum_separation_factor / 2);
                        disp_p2 *= (1. + minimum_separation_factor / 2);

                        std::pair<Real3, FaceID> pos_p1(p1.position(), fid1);
                        std::pair<Real3, FaceID> pos_p2(p2.position(), fid2);
                        pos_p1 = ecell4::polygon::travel(this->polygon(), pos_p1, disp_p1, 2);
                        pos_p2 = ecell4::polygon::travel(this->polygon(), pos_p2, disp_p2, 2);

                        p1.position() = pos_p1.first;
                        p2.position() = pos_p2.first;
                        this->update_particle(pid1, p1, pos_p1.second);
                        this->update_particle(pid2, p2, pos_p2.second);

                        boost::container::small_vector<
                            std::tuple<ParticleID, Particle, FaceID>, 2
                            > results(2);
                        results[0] = std::make_tuple(pid1, p1, pos_p1.second);
                        results[1] = std::make_tuple(pid2, p2, pos_p2.second);

                        STAT(stat_reaction_condition.add_count(PairSecondOrderFailed));
                        return results;
                    }
                }
                SGFRD_TRACE(tracer_.write("reaction is accepted."))
                this->update_particle(pid1, p_new, fid_new);

                this->remove_particle(pid2, fid2);
                last_reactions_.push_back(std::make_pair(rule,
                    make_binding_reaction_info(
                        tm, pid1, p1, pid2, p2, pid1, p_new)));

                // reaction succeed.
                STAT(stat_reaction_condition.add_count(PairSecondOrder));

                boost::container::small_vector<
                    std::tuple<ParticleID, Particle, FaceID>, 2> retval(1ul);
                retval[0] = std::make_tuple(pid1, p_new, fid_new);
                return retval;
            }
            default:
            {
                throw NotImplemented("SGFRD Pair Reaction: "
                    "more than two products from one reactant are not allowed");
            }
        }
    }

    Pair create_pair(const std::pair<ShellID, circle_type>& sh,
                     const ParticleID& pid1, const Particle& p1,
                     const ParticleID& pid2, const Particle& p2,
                     const Real3& ipv, const Real len_ipv)
    {
        SGFRD_SCOPE(ns, create_circular_pair_domain, tracer_);

        SGFRD_TRACE(tracer_.write("shell size = %1%", sh.second.size()))
        SGFRD_TRACE(tracer_.write("D1  = %1%, D2 = %2%", p1.D(), p2.D()))
        SGFRD_TRACE(tracer_.write("r1  = %1%, r2 = %2%", p1.radius(), p2.radius()))

        // ------------------ GF event (escape | ipv) --------------------------

        // the effective shell size is smaller than the actual shell size
        // by the radius of particles inside
        const Real shell_size = sh.second.size() - std::max(p1.radius(), p2.radius());

        const greens_functions::GreensFunction2DAbsSym
            gf_com(Pair::calc_D_com(p1.D(), p2.D()),
                   Pair::calc_R_com(shell_size, p1, p2));
        const Real t_com_escape = gf_com.drawTime(this->uniform_real());

        const Real k_tot = this->calc_k_tot(this->model_->query_reaction_rules(
                                           p1.species(), p2.species()));
        SGFRD_TRACE(tracer_.write("ipv length = %1%", len_ipv))
        SGFRD_TRACE(tracer_.write("total rate = %1%", k_tot))

        // TODO: [perf] if k_tot == 0.0 then we don't need GF.
        const greens_functions::GreensFunction2DRadAbs
            gf_ipv(Pair::calc_D_ipv(p1.D(), p2.D()),
                   k_tot, len_ipv, p1.radius() + p2.radius(),
                   Pair::calc_R_ipv(shell_size, p1, p2));
        const Real t_ipv_event = gf_ipv.drawTime(this->uniform_real());

        const std::pair<Real, Pair::EventKind> gf_event(
            (t_ipv_event < t_com_escape) ?
            std::make_pair(t_ipv_event,  Pair::IV_UNDETERMINED) :
            std::make_pair(t_com_escape, Pair::COM_ESCAPE));

        SGFRD_TRACE(tracer_.write("com escape time = %1%", t_com_escape))
        SGFRD_TRACE(tracer_.write("ipv event  time = %1%", t_ipv_event))
        SGFRD_TRACE(tracer_.write("GF  event  time = %1%", gf_event.first))

        // ------------------- single reaction event ---------------------------

        const Real t_single_reaction_1 = this->draw_reaction_time(
            this->calc_k_tot(this->model_->query_reaction_rules(p1.species())));
        const Real t_single_reaction_2 = this->draw_reaction_time(
            this->calc_k_tot(this->model_->query_reaction_rules(p2.species())));

        const std::pair<Real, Pair::EventKind> single_event(
            (t_single_reaction_1 < t_single_reaction_2) ?
            std::make_pair(t_single_reaction_1, Pair::SINGLE_REACTION_1) :
            std::make_pair(t_single_reaction_2, Pair::SINGLE_REACTION_2));

        SGFRD_TRACE(tracer_.write("particle 1 single reaction time = %1%", t_single_reaction_1))
        SGFRD_TRACE(tracer_.write("particle 2 single reaction time = %1%", t_single_reaction_2))
        SGFRD_TRACE(tracer_.write("single reaction time   = %1%", single_event.first))

        // --------------------------------------------------------------------

        if(gf_event.first < single_event.first)
        {
            SGFRD_TRACE(tracer_.write("gf event occurs first"))
            return Pair(gf_event.second, gf_event.first, this->time(), sh.first,
                        shell_size,
                        std::make_pair(pid1, p1),
                        std::make_pair(pid2, p2),
                        len_ipv, ipv, k_tot);
        }
        else
        {
            SGFRD_TRACE(tracer_.write("single reaction occurs first"))
            return Pair(single_event.second, single_event.first, this->time(),
                        sh.first,
                        shell_size,
                        std::make_pair(pid1, p1),
                        std::make_pair(pid2, p2),
                        len_ipv, ipv, k_tot);
        }
    }


//----------------------------------- multi ------------------------------------

    void fire_multi(Multi& dom, DomainID did)
    {
        SGFRD_SCOPE(us, fire_multi, tracer_);
        SGFRD_TRACE(tracer_.write("fire multi(%1%) for default dt(%2%)", did, dom.dt()));

        STAT(stat_multi_size.add_count(dom.particles().size());)

        volume_clearer vc(did, dom, *this, this->imm_sh_vis_applier);
        dom.step(vc);
        switch(dom.eventkind())
        {
            case Multi::NONE:
            {
                SGFRD_TRACE(tracer_.write("nothing occurs"))
                /* continuing multi domain: add this domain to scheduler */
                dom.begin_time() = this->time();
                this->add_event(dom);
                return;
            }
            case Multi::REACTION:
            {
                SGFRD_TRACE(tracer_.write("reaction occurs {");
                    for(auto rrec : dom.last_reactions())
                    {
                        tracer_.write("rule: %1% t: %2% reactant: { ",
                                      rrec.first.as_string(), rrec.second.t());
                        for(auto reactant : rrec.second.reactants())
                        {
                            tracer_.write("%1% ", reactant.first);
                        }
                        tracer_.write("} products: { ");
                        for(auto product: rrec.second.products())
                        {
                            tracer_.write("%1% ", product.first);
                        }
                        tracer_.write("}");
                    }
                    tracer_.write("}");)

                std::copy(dom.last_reactions().begin(), dom.last_reactions().end(),
                          std::back_inserter(this->last_reactions_));

                // record statistics
                for(const auto& rrec : dom.last_reactions())
                {
                    if(rrec.second.reactants().size() == 1)
                    {
                        STAT(stat_reaction_condition.add_count(MultiFirstOrder));
                    }
                    else if(rrec.second.reactants().size() == 2)
                    {
                        STAT(stat_reaction_condition.add_count(MultiSecondOrder));
                    }
                    else
                    {
                        throw std::runtime_error("unknown reaction record found");
                    }
                }

                ParticleID pid; Particle p; FaceID fid;
                for(const auto& pidpf : this->remove_multi(dom))
                {
                    std::tie(pid, p, fid) = pidpf;
                    this->add_event(this->create_tight_domain(
                        this->create_tight_shell(pid, p, fid), pid, p));
                }
                SGFRD_TRACE(tracer_.write("multi domain (id = %1%) removed.", did))
                return;
            }
            case Multi::ESCAPE:
            {
                SGFRD_TRACE(tracer_.write("escape occurs"))
                /* burst this domain! */
                ParticleID pid; Particle p; FaceID fid;
                for(const auto& pidpf : this->remove_multi(dom))
                {
                    std::tie(pid, p, fid) = pidpf;
                    this->add_event(this->create_tight_domain(
                        this->create_tight_shell(pid, p, fid), pid, p));
                }
                SGFRD_TRACE(tracer_.write("multi domain (id = %1%) removed.", did))
                return;
            }
            default:
            {
                SGFRD_TRACE(tracer_.write("Multi eventkind become invalid value!"))
                throw std::logic_error("never reach here");
            }
        }
    }

    // simply remove all the shells. not add a domains for each particles.
    // particles are updated at each step, so here nothing is needed to
    // update world.
    bursted_type remove_multi(const Multi& dom)
    {
        SGFRD_SCOPE(ns, remove_multi, tracer_);
        bursted_type results;
        Particle p; ParticleID pid;
        for(const auto& pidp : dom.particles())
        {
            std::tie(pid, p) = pidp;
            results.push_back(std::make_tuple(pid, p, this->get_face_id(pid)));
        }

        SGFRD_TRACE(tracer_.write("particles are collected"))

        for(const ShellID& sid : dom.shell_ids())
        {
            this->remove_shell(sid);
        }
        SGFRD_TRACE(tracer_.write("shells are removed"))
        return results;
    }

    // burst. step until(tm - dom.begin_time()).
    bursted_type burst_multi(Multi& dom, const Real tm)
    {
        SGFRD_SCOPE(ns, burst_multi, tracer_);

        auto did = get_domain_id(dom);
        volume_clearer vc(did, dom, *this, this->imm_sh_vis_applier);
        dom.step(vc, tm - dom.begin_time());

        SGFRD_TRACE(tracer_.write("multi domain steps with delta_t = %1%",
                    tm - dom.begin_time()));

        if(dom.eventkind() == Multi::REACTION)
        {
            SGFRD_TRACE(tracer_.write("reaction occured"));
            std::copy(dom.last_reactions().begin(), dom.last_reactions().end(),
                      std::back_inserter(this->last_reactions_));
        }
        return remove_multi(dom);
    }

//----------------------------------- birth ------------------------------------

    void fire_birth(const Birth& dom, DomainID did)
    {
        //TODO
        return;
    }


// -----------------------------------------------------------------------------

    // XXX: Note that the second element, distance to the domain, is
    //      - a distance to the Multi domains that ware not bursted or
    //      - a distance to the min-shell of particles that were bursted.
    std::vector<std::pair<DomainID, Real> >
    burst_and_shrink_non_multis(
            const ParticleID& pid, const Particle& p, const FaceID& fid,
            const std::vector<std::pair<DomainID, Real> >& intruders)
    {
        SGFRD_SCOPE(us, burst_and_shrink_non_multis, tracer_)

        const Real tm(this->time());
        std::vector<std::pair<DomainID, Real> > results;

        DomainID did; Real dist;
        for(const auto& didd : intruders)
        {
            SGFRD_SCOPE(ns, loop_for_intruders, tracer_)

            std::tie(did, dist) = didd;
            const auto& ev = get_event(did);
            if(ev->which_domain() == event_type::multi_domain)
            {
                SGFRD_TRACE(tracer_.write("domain %1% is multi", did))
                results.push_back(std::make_pair(did, dist));
                continue;
            }
            SGFRD_TRACE(tracer_.write("domain %1% is not a multi", did))

            DomainID did_; ParticleID pid_; Particle p_; FaceID fid_;
            for(const auto& pidpf : burst_event(std::make_pair(did, ev), tm))
            {
                std::tie(pid_, p_, fid_) = pidpf;

                SGFRD_TRACE(tracer_.write(
                    "add tight domain to bursted particle %1%", pid_))
                did_ = add_event(create_tight_domain(
                    create_tight_shell(pid_, p_, fid_), pid_, p_));

                // distance between the nearest one
                const auto dist = ecell4::polygon::distance(this->polygon(),
                        std::make_pair(p.position(), fid),
                        std::make_pair(p_.position(), fid_));

                // min shell radius of the nearest one
                const auto min_shell_rad =
                    calc_min_single_circular_shell_radius(p_) *
                    single_circular_shell_factor;

                // calculate distance as if the nearest one has a minimum shell
                results.emplace_back(did_, dist - min_shell_rad);
            }
            remove_event(did);
            SGFRD_TRACE(tracer_.write("domain %1% is bursted and shrinked", did))
        }

        std::sort(results.begin(), results.end(),
            ecell4::utils::pair_second_element_comparator<DomainID, Real>());

        SGFRD_TRACE(tracer_.write("results are sorted"))
        return results;
    }

    std::vector<std::pair<DomainID, Real> >
    burst_and_shrink_non_multis(const VertexID& vid,
        const std::vector<std::pair<DomainID, Real> >& intruders)
    {
        SGFRD_SCOPE(us, burst_and_shrink_non_multis_vertex, tracer_)
        const std::pair<Real3, VertexID> vpos = std::make_pair(
                polygon().position_at(vid), vid);

        const Real tm(this->time());
        std::vector<std::pair<DomainID, Real> > results;

        DomainID did; Real dist;
        for(const auto& didd : intruders)
        {
            SGFRD_SCOPE(ns, loop_for_intruders, tracer_)
            std::tie(did, dist) = didd;

            const auto& ev = get_event(did);
            if(ev->which_domain() == event_type::multi_domain)
            {
                SGFRD_TRACE(tracer_.write("domain %1% is multi", did))
                results.push_back(std::make_pair(did, dist));
                continue;
            }
            SGFRD_TRACE(tracer_.write("domain %1% is not a multi", did))

            DomainID did_; ParticleID pid_; Particle p_; FaceID fid_;
            for(const auto& pidpf : burst_event(std::make_pair(did, ev), tm))
            {
                std::tie(pid_, p_, fid_) = pidpf;
                SGFRD_TRACE(tracer_.write(
                    "add closely-fitted domain to bursted particle %1%", pid_))

                did_ = add_event(create_tight_domain(
                    create_tight_shell(pid_, p_, fid_), pid_, p_));

                const auto dist = ecell4::polygon::distance(this->polygon(),
                    vpos, std::make_pair(p_.position(), fid_));

                const auto min_shell_rad =
                    calc_min_single_circular_shell_radius(p_) *
                    single_circular_shell_factor;

                results.emplace_back(did_, dist - min_shell_rad);
            }
            remove_event(did);
            SGFRD_TRACE(tracer_.write("domain %1% is bursted and shrinked", did))
        }

        std::sort(results.begin(), results.end(),
            ecell4::utils::pair_second_element_comparator<DomainID, Real>());
        SGFRD_TRACE(tracer_.write("results are sorted"))
        return results;
    }


    // to clear volume. burst all the overlapping shells then add closely-fitted
    // shells to them. returns true if there are no overlapping particles.
    bool burst_and_shrink_overlaps(
            const Particle& p, const FaceID& fid, const DomainID& did);

    // form multi shell recursively
    DomainID form_multi(const ParticleID& pid, const Particle& p, const FaceID& fid,
                        const std::vector<std::pair<DomainID, Real> >& doms);

    // search intruder for multi, burst them and add them to multi if needed.
    void add_to_multi_recursive(Multi&);

    void merge_multi(Multi& from, Multi& to)
    {
        SGFRD_SCOPE(us, merge_multi, tracer_)

        const auto id_of_from = get_domain_id(from);
        remove_event(id_of_from);
        SGFRD_TRACE(tracer_.write("remove domain from(%1%)", id_of_from))

        // reset domain_id
        const domain_id_setter didset(this->get_domain_id(to));
        mut_sh_vis_applier(didset, from);
        SGFRD_TRACE(tracer_.write("set domain ID for shell in from(%1%)", id_of_from))


        // move particle
        ParticleID pid; Particle p;
        for(const auto& pidp : from.particles())
        {
            std::tie(pid, p) = pidp;
            const bool addp_result = to.add_particle(pid);
            assert(addp_result);
        }
        // move shell
        for(const ShellID& sid : from.shell_ids())
        {
            const bool adds_result = to.add_shell(sid);
            SGFRD_TRACE(tracer_.write("Shell(%1%).domain_id = %2%", sid,
                        boost::get<circular_shell_type>(this->get_shell(sid)).domain_id()))
            assert(adds_result);
        }
        to.determine_parameters();

        SGFRD_TRACE(tracer_.write("multi domain %1% is removed", id_of_from))
        return;
    }

    struct volume_clearer
    {
        volume_clearer(domain_id_type d, const Multi& dom, SGFRDSimulator& s,
                       immutable_shell_visitor_applier_type& imm)
            : sim(s), did(d), domain(dom), applier(imm)
        {}

        bool operator()(const Particle& p, const FaceID& fid)
        {
            SGFRD_SCOPE(us, volume_clearer, sim.tracer_)
            SGFRD_TRACE(sim.tracer_.write("ignoring nothing"))

            this->escaped_ = false;
            inside_checker is_inside(p.position(), p.radius(), fid, sim.polygon());
            if(applier(is_inside, domain))
            {
                SGFRD_TRACE(sim.tracer_.write("particle is inside of the shell"))
                return true;
            }
            SGFRD_TRACE(sim.tracer_.write("particle escaped. checking overlapping..."));

            const bool no_overlap = sim.burst_and_shrink_overlaps(p, fid, did);
            escaped_ = no_overlap;
            SGFRD_TRACE(sim.tracer_.write("no_overlap = %1%", no_overlap));
            return no_overlap;
        }
        bool operator()(const Particle& p, const FaceID& fid,
                        const ParticleID& ignore)
        {
            SGFRD_SCOPE(us, volume_clearer, sim.tracer_)
            SGFRD_TRACE(sim.tracer_.write("ignoring %1%", ignore))

            escaped_ = false;
            inside_checker is_inside(p.position(), p.radius(), fid, sim.polygon());
            if(applier(is_inside, domain))
            {
                SGFRD_TRACE(sim.tracer_.write("particle is inside of the shell"))
                return true;
            }
            SGFRD_TRACE(sim.tracer_.write("particle escaped. checking overlapping..."));

            const bool no_overlap = sim.burst_and_shrink_overlaps(p, fid, did);
            escaped_ = no_overlap;
            SGFRD_TRACE(sim.tracer_.write("no_overlap = %1%", no_overlap));
            return no_overlap;
        }
        bool operator()(const Particle& p, const FaceID& fid,
                        const ParticleID& ignore1, const ParticleID& ignore2)
        {
            SGFRD_SCOPE(us, volume_clearer, sim.tracer_)
            SGFRD_TRACE(sim.tracer_.write("ignoring %1% and %2%", ignore1, ignore2))

            escaped_ = false;
            inside_checker is_inside(p.position(), p.radius(), fid, sim.polygon());
            if(applier(is_inside, domain))
            {
                SGFRD_TRACE(sim.tracer_.write("particle is inside of the shell"))
                return true;
            }
            SGFRD_TRACE(sim.tracer_.write("particle escaped. checking overlapping..."));

            const bool no_overlap = sim.burst_and_shrink_overlaps(p, fid, did);
            escaped_ = no_overlap;
            SGFRD_TRACE(sim.tracer_.write("no_overlap = %1%", no_overlap));
            return no_overlap;
        }

        bool escaped() const {return escaped_;}

        tracer& access_tracer()
        {
            return sim.tracer_;
        }

      private:
        bool            escaped_;
        SGFRDSimulator& sim;
        domain_id_type  did;
        Multi const&    domain;
        immutable_shell_visitor_applier_type applier;
    };

//----------------------------------- event ------------------------------------

    //! make event from domain and push it into scheduler
    template<typename domainT>
    DomainID add_event(const domainT& dom)
    {
        SGFRD_SCOPE(ns, add_event, tracer_);
        SGFRD_TRACE(tracer_.write("the domain has dt = %1%, begin_time = %2%",
                    dom.dt(), dom.begin_time()))

        auto ev = boost::make_shared<event_type>(dom.begin_time() + dom.dt(), dom);
        const DomainID did = scheduler_.add(ev);

        SGFRD_TRACE(tracer_.write("event_time = %1%, domain ID = %2%", ev->time(), did))

        domain_id_setter didset(did);
        mut_sh_vis_applier(didset, dom);
        return did;
    }

    // assuming the event is already poped
    void fire_event(event_id_pair_type ev)
    {
        SGFRD_SCOPE(us, fire_event, tracer_);
        SGFRD_TRACE(tracer_.write("event %1% fired", ev.first))

        SGFRDEvent::domain_type& dom = ev.second->domain();

        switch(ev.second->which_domain())
        {
            case event_type::single_domain:
            {
                return this->fire_single(boost::get<Single>(dom), ev.first);
            }
            case event_type::pair_domain:
            {
                STAT(stat_fired_events.add_count(FirePair);)
                return this->fire_pair(boost::get<Pair>(dom), ev.first);
            }
            case event_type::multi_domain:
            {
                STAT(stat_fired_events.add_count(FireMulti);)
                return this->fire_multi(boost::get<Multi>(dom), ev.first);
            }
            case event_type::birth_domain:
            {
                STAT(stat_fired_events.add_count(FireBirth);)
                return this->fire_birth(boost::get<Birth>(dom), ev.first);
            }
            default:
            {
                throw std::runtime_error((boost::format(
                    "event::which_domain returns invalid value (%1%)") %
                    ev.second->which_domain()).str());
            }
        }
    }

    // assuming the event is already poped
    bursted_type burst_event(const event_id_pair_type& ev, Real tm)
    {
        SGFRD_SCOPE(us, burst_event, tracer_);
        SGFRDEvent::domain_type& dom = ev.second->domain();

        switch(ev.second->which_domain())
        {
            case event_type::single_domain:
            {
                return this->burst_single(boost::get<Single>(dom), tm);
            }
            case event_type::pair_domain:
            {
                return this->burst_pair(boost::get<Pair>(dom), tm);
            }
            case event_type::multi_domain:
            {
                return this->burst_multi(boost::get<Multi>(dom), tm);
            }
            default:
            {
                throw std::runtime_error((boost::format(
                    "event::which_domain returns invalid value (%1%)") %
                    ev.second->which_domain()).str());
            }
        }
    }

    ShellID create_tight_shell(
            const ParticleID& pid, const Particle& p, const FaceID fid)
    {
        SGFRD_SCOPE(us, create_tight_shell, tracer_);
        SGFRD_TRACE(tracer_.write("add close shell for %1% @ face %2%", pid, fid))
        SGFRD_TRACE(tracer_.write("species %1% has radius %2%",
                    p.species_serial(), p.radius()))

        const ShellID sid(shell_id_gen());
        circular_shell_type sh(circle_type(p.radius(), p.position(),
                               this->polygon().triangle_at(fid).normal()), fid);
        SGFRD_TRACE(tracer_.write("shell has size == %1%", sh.size()))

        shell_container_.check_add_shell(sid, sh, fid, "create tight shell");
        SGFRD_TRACE(tracer_.write("new shell id is %1%", sid))
        return sid;
    }
    Single create_tight_domain(
            const ShellID& sid, const ParticleID& pid, const Particle& p)
    {
        SGFRD_SCOPE(us, create_tight_domain, tracer_);
        SGFRD_TRACE(tracer_.write("for particle %1%, shell %2%", pid, sid))
        return Single(Single::ESCAPE, 0., this->time(), sid, std::make_pair(pid, p));
    }

    std::pair<ShellID, circle_type>
    create_minimum_single_shell(
            const ParticleID& pid, const Particle& p, const FaceID fid)
    {
        // this function is used to create shells for multi.
        // multi domain allows its shells to overlap each other.
        return create_single_circular_shell(std::make_pair(p.position(), fid),
            calc_min_single_circular_shell_radius(p), false);
    }

    Multi create_empty_multi()
    {
        return Multi(*this, *(this->world_), this->time(),
                     this->bd_dt_factor_, this->reaction_length_);
    }

    //! make domain and call add_event
    DomainID create_event(const ParticleID&, const Particle&, const FaceID);

    expected<DomainID, std::vector<std::pair<DomainID, Real> > >
    form_single_circular_event(
            const ParticleID&, const Particle&, const FaceID, const Real);

    expected<DomainID, std::vector<std::pair<DomainID, Real> > >
    form_single_conical_event(
            const ParticleID&, const Particle&, const FaceID);

    boost::optional<DomainID>
    form_pair(const ParticleID& pid, const Particle& p, const FaceID& fid,
              const std::vector<std::pair<DomainID, Real> >& intruders);


    std::vector<std::pair<VertexID, Real> >
    get_intrusive_vertices(const std::pair<Real3, FaceID>& pos,
                           const Real radius) const
    {
        return this->world_->list_vertices_within_radius(pos, radius);
    }

    // the second value is distance between
    // the point `pos` and the surface of a domain
    std::vector<std::pair<DomainID, Real>>
    get_intrusive_domains(const std::pair<Real3, FaceID>& pos,
                          const Real radius) const
    {
        SGFRD_SCOPE(us, get_intrusive_domains_position, tracer_);

        const auto shells(shell_container_.list_shells_within_radius(pos, radius));

        SGFRD_TRACE(tracer_.write("collected %1% shells in radius %2%",
                    shells.size(), radius));

        std::vector<std::pair<DomainID, Real>> domains;
        domains.reserve(shells.size());
        for(const auto& shid_sh : shells)
        {
            const auto& shell_id_pair = shid_sh.first;
            const auto  dist          = shid_sh.second;

            SGFRD_TRACE(tracer_.write("shell %1% = {%2%} is at %3% distant",
                        shell_id_pair.first, shell_id_pair.second, dist));

            const DomainID did = boost::apply_visitor(
                    domain_id_getter(), shell_id_pair.second);

            SGFRD_TRACE(tracer_.write("shell %1% is related to domain %2%",
                        shell_id_pair.first, did));

            // filter by domainID to be unique
            auto found = std::find_if(domains.begin(), domains.end(),
                [did](const std::pair<DomainID, Real>& did_dist) noexcept -> bool {
                    return did_dist.first == did;
                });
            if(found == domains.end())
            {
                SGFRD_TRACE(tracer_.write("domain %1% is assigned to retval", did));
                domains.emplace_back(did, dist);
            }
            else
            {
                // choose the nearer one from the domain that contains multiple
                // shells.
                found->second = std::min(found->second, dist);
                SGFRD_TRACE(tracer_.write("domain %1% is already assigned", did));
            }
        }
        for(std::size_t i=0; i < domains.size(); ++i)
        {
            SGFRD_TRACE(tracer_.write("%1%, ", domains.at(i).first));
        }

        std::sort(domains.begin(), domains.end(),
                  utils::pair_second_element_comparator<DomainID, Real>());
        return domains;
    }

    std::vector<std::pair<DomainID, Real> >
    get_intrusive_domains(const VertexID& vid, const Real radius) const
    {
        SGFRD_SCOPE(us, get_intrusive_domains_vid, tracer_);

        const std::pair<Real3, VertexID> vpos(polygon().position_at(vid), vid);
        const std::vector<std::pair<std::pair<ShellID, shell_type>, Real>
            > shells(shell_container_.list_shells_within_radius(vpos, radius));

        SGFRD_TRACE(tracer_.write("collected %1% shells in radius %2%",
                    shells.size(), radius));

        std::vector<std::pair<DomainID, Real> > domains;
        domains.reserve(shells.size());

        for(std::vector<std::pair<std::pair<ShellID, shell_type>, Real>
                >::const_iterator iter(shells.begin()), iend(shells.end());
                iter != iend; ++iter)
        {
            std::pair<ShellID, shell_type> shell_id_pair;
            Real dist;
            std::tie(shell_id_pair, dist) = *iter;
            SGFRD_TRACE(tracer_.write("shell %1% = {%2%} is at %3% distant",
                        shell_id_pair.first, shell_id_pair.second, dist));

            const DomainID did = boost::apply_visitor(
                    domain_id_getter(), shell_id_pair.second);

            SGFRD_TRACE(tracer_.write("shell %1% is related to domain %2%",
                        shell_id_pair.first, did));

            std::vector<std::pair<DomainID, Real> >::iterator found =
                std::find_if(domains.begin(), domains.end(),
                    ecell4::utils::pair_first_element_unary_predicator<
                    DomainID, Real>(did));
            if(found == domains.end())
            {
                SGFRD_TRACE(tracer_.write("domain %1% is assigned to retval", did));
                domains.push_back(std::make_pair(did, dist));
            }
            else
            {
                found->second = std::min(found->second, dist);
                SGFRD_TRACE(tracer_.write("domain %1% is already assigned", did));
                for(std::size_t i=0; i < domains.size(); ++i)
                {
                    SGFRD_TRACE(tracer_.write("%1%, ", domains.at(i).first));
                }
            }
        }
        std::sort(domains.begin(), domains.end(),
                  utils::pair_second_element_comparator<DomainID, Real>());

        return domains;
    }

    //! just a geometric restriction
    Real get_max_circle_size(const std::pair<Real3, FaceID>& pos) const
    {
        SGFRD_SCOPE(us, get_max_circle_size, tracer_);
        SGFRD_TRACE(tracer_.write("position = %1%, %2%", pos.first, pos.second));

        Real lensq = std::numeric_limits<Real>::max();
        for(const auto& barrier : world_->barrier_at(pos.second))
        {
            const Real dist2 = ecell4::sgfrd::distance_sq(pos.first, barrier);
            if(dist2 < lensq)
            {
                lensq = dist2;
            }
        }
        SGFRD_TRACE(tracer_.write("distance to barrier = %1%", std::sqrt(lensq)));
        return std::sqrt(lensq);
    }
    Real get_max_cone_size(const VertexID& vid) const
    {
        Real min_len = std::numeric_limits<Real>::max();
        for(const auto eid : this->polygon().outgoing_edges(vid))
        {
            min_len = std::min(min_len, this->polygon().length_of(eid));
        }
        return min_len * 0.5;
    }

    Real calc_min_single_circular_shell_radius(const Particle& p) const
    {
        return p.radius() * single_circular_shell_factor;
    }
    Real calc_min_single_conical_shell_radius(const Particle& p) const
    {
        return p.radius() * single_conical_surface_shell_factor;
    }

    Real calc_k_tot(const std::vector<ReactionRule>& rules) const
    {
        SGFRD_SCOPE(ns, calc_k_tot, tracer_)
        Real k_tot = 0;
        for(const auto& rule : rules)
        {
            SGFRD_TRACE(tracer_.write("reaction rule found: k = %1%", rule.k()));
            k_tot += rule.k();
        }
        return k_tot;
    }

    Real draw_reaction_time(const Real k_tot)
    {
        SGFRD_SCOPE(ns, draw_reaction_time, tracer_)
        SGFRD_TRACE(tracer_.write("for ktot = %1%", k_tot));

        if(k_tot <= 0){return std::numeric_limits<Real>::infinity();}
        if(k_tot == std::numeric_limits<Real>::infinity()){return 0;}

        const Real rnd = this->uniform_real();
        SGFRD_TRACE(tracer_.write("rnd = %1%", rnd));
        if(rnd <= 0){return std::numeric_limits<Real>::infinity();}

        SGFRD_TRACE(tracer_.write("reaction occurs in a finite time"));

        return (1. / k_tot) * (-std::log(rnd));
    }

    ReactionRule const&
    determine_reaction_rule(const std::vector<ReactionRule>& rules)
    {
        if(rules.size() == 1){return rules.front();}
        if(rules.empty())
        {
            throw std::invalid_argument("no reaction rule exists");
        }

        const Real ktot = calc_k_tot(rules);
        const Real thrs = ktot * this->uniform_real();
        Real a = 0.0;
        for(const auto& rule : rules)
        {
            a += rule.k();
            if(a > thrs) return rule;
        }
        throw std::logic_error("reaction cannot detemined");
    }

    // for BD part...
    Real3 random_circular_uniform(const Real r, const FaceID& fid)
    {
        SGFRD_SCOPE(ns, sgfrd_random_circular_uniform, tracer_);

        const Real   theta = this->rng_.uniform(0., 2 * M_PI);
        SGFRD_TRACE(tracer_.write("theta = %1%", theta));

        const Real3  rnd(r * std::cos(theta), r * std::sin(theta), 0.);
        SGFRD_TRACE(tracer_.write("random displacement = %1%", rnd));

        const Real3& normal = this->polygon().triangle_at(fid).normal();
        SGFRD_TRACE(tracer_.write("normal vector on face(%1%) = %2%", fid, normal));

        const Real   tilt   = calc_angle(Real3(0, 0, 1), normal);
        SGFRD_TRACE(tracer_.write("tilt angle = %1%", tilt));

        if     (std::abs(tilt)        < 1e-10) {return rnd;}
        else if(std::abs(tilt - M_PI) < 1e-10) {return rnd * (-1.0);}
        else if(std::abs(tilt + M_PI) < 1e-10) {return rnd * (-1.0);}

        SGFRD_TRACE(tracer_.write("tilt angle is not 0 or +-pi. rotating..."));

        const Real3 axis = cross_product(Real3(0., 0., 1.), normal);
        return rotate(tilt, axis * (1. / length(axis)), rnd);
    }

    static Real calc_modest_shell_size(
            const Particle& p, const Particle& nearest, const Real dist)
    {
        // here, `dist` means the distance between shell surfaces.
        // so the following relationship should be satisfied.
        // >> p.radius() + nearest.radius() + dist ==
        // >> distance_on_polygon(p.position(), nearest.position())

        // However, this function does not take the face information,
        // we cannot calculate the distance between p and nearest.
        assert(0.0 <= dist);

        const Real D1     = p.D();
        const Real D2     = nearest.D();
        assert(!(D1 == 0.0 && D2 == 0.0));
        const Real r1     = p.radius();
        const Real sqrtD1 = std::sqrt(D1);
        const Real sqrtD2 = std::sqrt(D2);
        return sqrtD1 / (sqrtD1 + sqrtD2) * dist + r1;
    }

  private:

    static const Real single_circular_shell_factor;
    static const Real single_circular_shell_mergin;
    static const Real single_conical_surface_shell_factor;
    static const Real single_conical_surface_shell_mergin;
    static const Real minimum_separation_factor;

  private:

    // from SimulatorBase
    // boost::shared_ptr<model_type> model_;
    // boost::shared_ptr<world_type> world_;
    // Integer num_steps_;
    bool is_dirty_;
    Real dt_;
    Real bd_dt_factor_;
    Real reaction_length_;
    ecell4::RandomNumberGenerator&       rng_;
    scheduler_type                       scheduler_;
    shell_id_generator_type              shell_id_gen;
    shell_container_type                 shell_container_;
    mutable_shell_visitor_applier_type   mut_sh_vis_applier;
    immutable_shell_visitor_applier_type imm_sh_vis_applier;
    std::vector<std::pair<reaction_rule_type, reaction_info_type> > last_reactions_;
    mutable tracer tracer_;

  public:
    STAT(mutable statistics<ReactionKind> stat_reaction_condition;)
    STAT(mutable statistics<MultiReason> stat_multi_reason;)
    STAT(mutable statistics<EventFired>  stat_fired_events;)
    STAT(mutable statistics<std::size_t> stat_multi_size;)
};


} // sgfrd
} // ecell4
#endif // ECELL4_SGFRD_SIMULATOR
