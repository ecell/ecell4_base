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
#include <ecell4/sgfrd/Informations.hpp>
#include <ecell4/sgfrd/SGFRDEvent.hpp>
#include <ecell4/sgfrd/expected.hpp>

#include <ecell4/sgfrd/tracer.hpp>

#include <boost/make_shared.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/container/static_vector.hpp>
#include <boost/container/small_vector.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>

#include <iostream>

#define SGFRD_LOG(x, y) /**/

namespace ecell4
{
namespace sgfrd
{

class SGFRDSimulator :
    public ecell4::SimulatorBase<ecell4::Model, SGFRDWorld>
{
  public:
    typedef SGFRDSimulator self_type;

    // polygon
    typedef polygon_traits polygon_traits_type;
    typedef ecell4::Polygon<polygon_traits_type>     polygon_type;
    typedef polygon_type::triangle_type     triangle_type;
    typedef polygon_type::face_id_type      face_id_type;
    typedef polygon_type::edge_id_type      edge_id_type;
    typedef polygon_type::vertex_id_type    vertex_id_type;
    typedef polygon_type::face_descripter   face_descripter;
    typedef polygon_type::edge_descripter   edge_descripter;
    typedef polygon_type::vertex_descripter vertex_descripter;
    typedef polygon_type::local_index_type  local_index_type;
    typedef polygon_type::barycentric_type  barycentric_type;
    typedef face_id_type   FaceID;// just for looks same as ParticleID
    typedef edge_id_type   EdgeID;
    typedef vertex_id_type VertexID;

    // Event & Domain
    typedef SGFRDEvent                      event_type;
    typedef SGFRDEvent::domain_type         domain_type;
    typedef EventID                         event_id_type;
    typedef DomainID                        domain_id_type;
    typedef SGFRDEventScheduler             scheduler_type;
    typedef SGFRDEventScheduler::value_type event_id_pair_type;

    // Simulator
    typedef ecell4::SimulatorBase<ecell4::Model, SGFRDWorld> base_type;
    typedef base_type::world_type world_type;
    typedef base_type::model_type model_type;
    typedef std::pair<ParticleID, Particle> particle_id_pair_type;
    typedef boost::tuple<ParticleID, Particle, FaceID> pid_p_fid_tuple_type;
    // length of this type may be tuned
    typedef boost::container::small_vector<pid_p_fid_tuple_type, 5> bursted_type;

    // ShellContainer
    typedef ecell4::SerialIDGenerator<ShellID> shell_id_generator_type;
    typedef ShellContainer<polygon_traits_type> shell_container_type;
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
    typedef MoleculeInfo                   molecule_info_type;
    typedef ReactionInfo                   reaction_info_type;
    typedef std::pair<reaction_rule_type, reaction_info_type> reaction_log_type;
    typedef std::vector<reaction_log_type> reaction_archive_type;

  public:

    SGFRDSimulator(const boost::shared_ptr<world_type>& world,
                   const boost::shared_ptr<model_type>& model,
                   Real bd_dt_factor = 1e-5, Real reaction_length = 1e-3,
                   const std::string& trace_fname = "sgfrd_trace.log")
        : base_type(model, world), dt_(0),
          bd_dt_factor_(bd_dt_factor), reaction_length_(reaction_length),
          rng_(*(world->rng())), shell_container_(*(world->polygon())),
          mut_sh_vis_applier(shell_container_), imm_sh_vis_applier(shell_container_),
          tracer_(trace_fname)
    {}

    SGFRDSimulator(boost::shared_ptr<world_type> world,
                   Real bd_dt_factor = 1e-5, Real reaction_length = 1e-3,
                   const std::string& trace_fname = "sgfrd_trace.log")
        : base_type(world), dt_(0),
          bd_dt_factor_(bd_dt_factor), reaction_length_(reaction_length),
          rng_(*(world->rng())), shell_container_(*(world->polygon())),
          mut_sh_vis_applier(shell_container_), imm_sh_vis_applier(shell_container_),
          tracer_(trace_fname)
    {}

    ~SGFRDSimulator(){}

    void initialize()
    {
        ParticleID pid; Particle p;
        BOOST_FOREACH(boost::tie(pid, p), this->world_->list_particles())
        {
            add_event(create_closely_fitted_domain(create_closely_fitted_shell(
                      pid, p, this->get_face_id(pid)), pid, p));
        }
        return ;
    }
    void finalize()
    {
        const Real tm(this->time());
        while(this->scheduler_.size() != 0)
        {
            this->burst_event(this->scheduler_.pop(), tm);
        }
        return ;
    }
    void step()
    {
        SGFRD_SCOPE(us, step, tracer_);

        this->set_time(this->scheduler_.next_time());
        SGFRD_TRACE(tracer_.write("now t = %1%", this->time()))

        // fire event executes `create_event` inside.
        this->fire_event(this->scheduler_.pop());
        SGFRD_TRACE(tracer_.write("now %1% shells exists", shell_container_.num_shells()))
        SGFRD_TRACE(tracer_.write("now %1% events exists", scheduler_.size()))

        //XXX
        assert(this->diagnosis());

        return;
    }
    bool step(const Real& upto)
    {
        this->step();
        return this->time() < upto;
    }

    Real dt() const {return dt_;}
    Real reaction_length() const {return reaction_length_;}

    bool check_reaction() const {return last_reactions_.size() > 0;}
    std::vector<std::pair<ReactionRule, reaction_info_type> > const&
    last_reactions() const {return last_reactions_;}

    bool diagnosis() const;

  private:

    // simple wrappers to call member's member-method (e.g. world_->t()) {{{
    Real uniform_real(){return this->rng_.random();}

    world_type   const& world()   const {return *(this->world_);}
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
    void set_time(const Real t) {return this->world_->set_t(t);}

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

    boost::tuple<ParticleID, Particle, FaceID>
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

    boost::tuple<ParticleID, Particle, FaceID>
    propagate_single_circular(
        const circular_shell_type& sh, const Single& dom, const Real tm);
    boost::tuple<ParticleID, Particle, FaceID>
    propagate_single_conical(
        const conical_surface_shell_type& sh, const Single& dom, const Real tm);

    boost::tuple<ParticleID, Particle, FaceID>
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

    boost::tuple<ParticleID, Particle, FaceID>
    escape_single_circular(const circular_shell_type& sh, const Single& dom);
    boost::tuple<ParticleID, Particle, FaceID>
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
            const std::pair<Real3, FaceID>& pos, const Real size)
    {
        SGFRD_SCOPE(ns, create_single_circular_shell, tracer_);

        const ShellID id(shell_id_gen());
        const circle_type shape(size, pos.first,
                                polygon().triangle_at(pos.second).normal());
        shell_container_.add_shell(id, circular_shell_type(shape, pos.second),
                                   pos.second);
        SGFRD_TRACE(tracer_.write("the shell id is %1%", id))
        return std::make_pair(id, shape);
    }
    std::pair<ShellID, conical_surface_type>
    create_single_conical_surface_shell(
            const vertex_id_type& vid, const Real size)
    {
        SGFRD_SCOPE(ns, create_single_conical_shell, tracer_);

        const ShellID id(shell_id_gen());
        const conical_surface_type shape(polygon().vertex_at(vid).position,
                                         polygon().apex_angle(vid), size);
        shell_container_.add_shell(
                id, conical_surface_shell_type(shape, vid), vid);
        SGFRD_TRACE(tracer_.write("the shell id is %1%", id));

        return std::make_pair(id, shape);
    }

    Single create_single(const std::pair<ShellID, circle_type>& sh,
                         const ParticleID& pid, const Particle& p)
    {
        SGFRD_SCOPE(ns, create_single_circular_domain, tracer_);

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

        SGFRD_TRACE(tracer_.write("shell size = %1%", sh.second.size()))
        SGFRD_TRACE(tracer_.write("D   = %1%", p.D()))
        SGFRD_TRACE(tracer_.write("r0  = %1%", length(p.position() - sh.second.apex())))
        SGFRD_TRACE(tracer_.write("a   = %1%", sh.second.size() - p.radius()))
        SGFRD_TRACE(tracer_.write("phi = %1%", sh.second.apex_angle()))

        const greens_functions::GreensFunction2DRefWedgeAbs
            gf(/* D   = */ p.D(),
               /* r0  = */ length(p.position() - sh.second.apex()),
               /* a   = */ sh.second.size() - p.radius(),
               /* phi = */ sh.second.apex_angle());
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
            else
            {
                dom.eventkind() = Pair::IV_REACTION;
                SGFRD_TRACE(tracer_.write("pair event kind = IV_REACTION"));
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
                boost::array<boost::tuple<ParticleID, Particle, FaceID>, 2>
                    propagated = this->propagate_pair(
                        this->get_shell(sid), dom, this->time());
                SGFRD_TRACE(tracer_.write("particles are propagated"));

                this->remove_shell(sid);
                SGFRD_TRACE(tracer_.write("shell %1% removed", sid));

                boost::array<ShellID,    2> sids;
                boost::array<DomainID,   2> dids;
                boost::array<Single,     2> doms;
                boost::array<ParticleID, 2> pids;
                boost::array<Particle,   2> ps;
                boost::array<FaceID,     2> fids;

                // add tight-domain for them to detect overlap
                for(std::size_t i=0; i<2; ++i)
                {
                    boost::tie(pids[i], ps[i], fids[i]) = propagated[i];
                    SGFRD_TRACE(tracer_.write("adding tight domain > %1%", pids[i]))
                    sids[i] = create_closely_fitted_shell(pids[i], ps[i], fids[i]);
                    doms[i] = create_closely_fitted_domain(sids[i], pids[i], ps[i]);
                    dids[i] = add_event(doms[i]);
                }
                SGFRD_TRACE(tracer_.write("tight-domains assigned"));

                // after that, attempt single reaction
                const std::size_t ridx = *reactant_index;
                SGFRD_TRACE(tracer_.write("reactant_index = %1%", ridx));

                BOOST_AUTO(results, this->attempt_reaction_single(
                           this->get_shell(sids[ridx]), dids[ridx], doms[ridx],
                           pids[ridx], ps[ridx], fids[ridx]));

                this->remove_shell(sids[ridx]);
                SGFRD_TRACE(tracer_.write("shell %1% removed", sids[ridx]));
                this->remove_event(dids[ridx]);
                SGFRD_TRACE(tracer_.write("event %1% removed", dids[ridx]));

                SGFRD_TRACE(tracer_.write("reaction attempted"));

                // add domain to each reactant
                ParticleID pid; Particle p; FaceID fid;
                BOOST_FOREACH(boost::tie(pid, p, fid), results)
                {
                    SGFRD_TRACE(tracer_.write("adding next event for %1%", pid))
                    add_event(create_closely_fitted_domain(
                              create_closely_fitted_shell(pid, p, fid), pid, p));
                }
                return;
            }
            case Pair::COM_ESCAPE:
            case Pair::IV_ESCAPE:
            {
                SGFRD_SCOPE(ns, case_COM_OR_IV_ESCAPE, tracer_);
                boost::array<boost::tuple<ParticleID, Particle, FaceID>, 2>
                    escaped = this->escape_pair(
                        this->get_shell(sid), dom, this->time());
                SGFRD_TRACE(tracer_.write("particles escaped"));

                this->remove_shell(sid);
                SGFRD_TRACE(tracer_.write("shell %1% removed", sid));

                ParticleID pid; Particle p; FaceID fid;
                BOOST_FOREACH(boost::tie(pid, p, fid), escaped)
                {
                    SGFRD_TRACE(tracer_.write("adding next event for %1%", pid))
                    add_event(create_closely_fitted_domain(
                                create_closely_fitted_shell(pid, p, fid), pid, p));
                }
                return;
            }
            case Pair::IV_REACTION:
            {
                SGFRD_SCOPE(ns, case_IV_REACTION, tracer_);
                boost::container::small_vector<
                    boost::tuple<ParticleID, Particle, FaceID>, 2>
                        products = this->attempt_pair_reaction(
                                this->get_shell(sid), dom, this->time());
                SGFRD_TRACE(tracer_.write("iv reaction attempted"));
                this->remove_shell(sid);
                SGFRD_TRACE(tracer_.write("shell %1% removed", sid));

                if(!products.empty())
                {
                    ParticleID pid; Particle p; FaceID fid;
                    boost::tie(pid, p, fid) = products.front();

                    SGFRD_TRACE(tracer_.write(
                                "adding next event for particle %1%", pid))

                    add_event(create_closely_fitted_domain(
                                create_closely_fitted_shell(pid, p, fid), pid, p));
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
        boost::array<boost::tuple<ParticleID, Particle, FaceID>, 2>
            propagated(this->propagate_pair(this->get_shell(sid), dom, tm));
        results.push_back(propagated[0]);
        results.push_back(propagated[1]);

        SGFRD_TRACE(tracer_.write("particle %1% and %2% propagated",
            boost::get<0>(propagated[0]), boost::get<0>(propagated[1])));

        this->remove_shell(sid);
        SGFRD_TRACE(tracer_.write("shell(%1%) removed", sid));
        return results;
    }

    boost::array<boost::tuple<ParticleID, Particle, FaceID>, 2>
    propagate_pair(const shell_type& sh, const Pair& dom, const Real tm)
    {
        SGFRD_SCOPE(us, propagate_pair, tracer_);
        switch(sh.which())
        {
            case shell_container_type::circular_shell:
            {
                return propagate_pair_circular(
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

    boost::array<boost::tuple<ParticleID, Particle, FaceID>, 2>
    propagate_pair_circular(const circular_shell_type& sh,
            const Pair& dom, const Real tm)
    {
        SGFRD_SCOPE(us, propagate_pair_circular, tracer_);

        // calculate displacements
        const Real dt = tm - dom.begin_time();
        const greens_functions::GreensFunction2DRadAbs
            gf_ipv(dom.D_ipv(), dom.kf(), dom.r0(), dom.sigma(), dom.R_ipv());
        const Real l_ipv     = gf_ipv.drawR(this->uniform_real(), dt);
        const Real theta_ipv = gf_ipv.drawTheta(this->uniform_real(), l_ipv, dt);

        const greens_functions::GreensFunction2DAbsSym
            gf_com(dom.D_com(), dom.R_com());
        const Real l_com     = gf_com.drawR(this->uniform_real(), dt);
        const Real theta_com = this->uniform_real() *
                               boost::math::constants::two_pi<Real>();

        const FaceID       sh_fid = sh.structure_id();
        const triangle_type&    f = this->polygon().triangle_at(sh_fid);
        const Real3 direction_com = rotate(theta_com, f.normal(), f.represent());

        const Real3 disp_com = direction_com * (l_com / length(direction_com));
        const Real3 disp_ipv = rotate(theta_ipv, f.normal(), dom.ipv());

        // update position
        // XXX to treat polygon surface structure, calculate displacement for
        // each partile without polygon information first
        // and then apply it to each particle with structure

        Particle         p1   = dom.particle_at(0);
        Particle         p2   = dom.particle_at(1);
        const ParticleID pid1 = dom.particle_id_at(0);
        const ParticleID pid2 = dom.particle_id_at(1);

        // ipv is a vector from p1 to p2
        const Real  ratio_p1     = -p1.D() / (p1.D() + p2.D());
        const Real  ratio_p2     =  p2.D() / (p1.D() + p2.D());
        const Real3 disp_ipv_p1  = disp_ipv  * ratio_p1;
        const Real3 disp_ipv_p2  = disp_ipv  * ratio_p2;
        const Real3 disp_ipv0_p1 = dom.ipv() * ratio_p1;
        const Real3 disp_ipv0_p2 = dom.ipv() * ratio_p2;

        Real3 disp_p1 = disp_com - disp_ipv0_p1 + disp_ipv_p1;
        Real3 disp_p2 = disp_com - disp_ipv0_p2 + disp_ipv_p2;
        std::pair<Real3, FaceID> pos_p1(p1.position(), this->get_face_id(pid1));
        std::pair<Real3, FaceID> pos_p2(p2.position(), this->get_face_id(pid2));

        if(0 == ecell4::polygon::travel(this->polygon(), pos_p1, disp_p1, 2))
        {
            SGFRD_TRACE(tracer_.write("escape_com_circular_pair "
                        "p1 moving on face: precision lost"))
        }
        if(0 == ecell4::polygon::travel(this->polygon(), pos_p2, disp_p2, 2))
        {
            SGFRD_TRACE(tracer_.write("escape_com_circular_pair "
                        "p2 moving on face: precision lost"))
        }

        p1.position() = pos_p1.first;
        p2.position() = pos_p2.first;

        this->update_particle(pid1, p1, pos_p1.second);
        this->update_particle(pid2, p2, pos_p2.second);

        boost::array<boost::tuple<ParticleID, Particle, FaceID>, 2> results;
        results[0] = boost::make_tuple(pid1, p1, pos_p1.second);
        results[1] = boost::make_tuple(pid2, p2, pos_p2.second);

        return results;
    }

    boost::array<boost::tuple<ParticleID, Particle, FaceID>, 2>
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
    boost::array<boost::tuple<ParticleID, Particle, FaceID>, 2>
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
    boost::array<boost::tuple<ParticleID, Particle, FaceID>, 2>
    escape_com_circular_pair(
            const circular_shell_type& sh, const Pair& dom, const Real tm)
    {
        SGFRD_SCOPE(us, escape_com_pair, tracer_);

        // calculate displacements
        const Real dt = tm - dom.begin_time();
        const greens_functions::GreensFunction2DRadAbs
            gf_ipv(dom.D_ipv(), dom.kf(), dom.r0(), dom.sigma(), dom.R_ipv());
        const Real l_ipv     = gf_ipv.drawR    (this->uniform_real(), dt);
        const Real theta_ipv = gf_ipv.drawTheta(this->uniform_real(), l_ipv, dt);
        const Real l_com     = dom.R_com();
        const Real theta_com = this->uniform_real() *
                               boost::math::constants::two_pi<Real>();

        const FaceID       sh_fid = sh.structure_id();
        const triangle_type&    f = this->polygon().triangle_at(sh_fid);
        const Real3 direction_com = rotate(theta_com, f.normal(), f.represent());

        const Real3 disp_com = direction_com * (l_com / length(direction_com));
        const Real3 disp_ipv = rotate(theta_ipv, f.normal(), dom.ipv());

        // update position
        // XXX to treat polygon surface structure, calculate displacement for
        // each partile without polygon information first
        // and then apply it to each particle with structure

        Particle         p1   = dom.particle_at(0);
        Particle         p2   = dom.particle_at(1);
        const ParticleID pid1 = dom.particle_id_at(0);
        const ParticleID pid2 = dom.particle_id_at(1);

        // ipv is a vector from p1 to p2
        const Real  ratio_p1     = -p1.D() / (p1.D() + p2.D());
        const Real  ratio_p2     =  p2.D() / (p1.D() + p2.D());
        const Real3 disp_ipv_p1  = disp_ipv  * ratio_p1;
        const Real3 disp_ipv_p2  = disp_ipv  * ratio_p2;
        const Real3 disp_ipv0_p1 = dom.ipv() * ratio_p1;
        const Real3 disp_ipv0_p2 = dom.ipv() * ratio_p2;

        Real3 disp_p1 = disp_com - disp_ipv0_p1 + disp_ipv_p1;
        Real3 disp_p2 = disp_com - disp_ipv0_p2 + disp_ipv_p2;
        std::pair<Real3, FaceID> pos_p1(p1.position(), this->get_face_id(pid1));
        std::pair<Real3, FaceID> pos_p2(p2.position(), this->get_face_id(pid2));

        if(0 == ecell4::polygon::travel(this->polygon(), pos_p1, disp_p1, 2))
        {
            SGFRD_TRACE(tracer_.write("escape_com_circular_pair "
                        "p1 moving on face: precision lost"))
        }
        if(0 == ecell4::polygon::travel(this->polygon(), pos_p2, disp_p2, 2))
        {
            SGFRD_TRACE(tracer_.write("escape_com_circular_pair "
                        "p1 moving on face: precision lost"))
        }

        p1.position() = pos_p1.first;
        p2.position() = pos_p2.first;

        this->update_particle(pid1, p1, pos_p1.second);
        this->update_particle(pid2, p2, pos_p2.second);

        boost::array<boost::tuple<ParticleID, Particle, FaceID>, 2> results;
        results[0] = boost::make_tuple(pid1, p1, pos_p1.second);
        results[1] = boost::make_tuple(pid2, p2, pos_p2.second);

        return results;
    }

    boost::array<boost::tuple<ParticleID, Particle, FaceID>, 2>
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
    boost::array<boost::tuple<ParticleID, Particle, FaceID>, 2>
    escape_ipv_circular_pair(
            const circular_shell_type& sh, const Pair& dom, const Real tm)
    {
        SGFRD_SCOPE(us, escape_com_pair, tracer_);

        // calculate displacements
        const Real dt = tm - dom.begin_time();
        const greens_functions::GreensFunction2DRadAbs
            gf_ipv(dom.D_ipv(), dom.kf(), dom.r0(), dom.sigma(), dom.R_ipv());
        const Real l_ipv     = dom.R_ipv();
        const Real theta_ipv = gf_ipv.drawTheta(this->uniform_real(), l_ipv, dt);

        const greens_functions::GreensFunction2DAbsSym
            gf_com(dom.D_com(), dom.R_com());
        const Real l_com     = gf_com.drawR(this->uniform_real(), dt);
        const Real theta_com = this->uniform_real() *
                               boost::math::constants::two_pi<Real>();

        const FaceID       sh_fid = sh.structure_id();
        const triangle_type&    f = this->polygon().triangle_at(sh_fid);
        const Real3 direction_com = rotate(theta_com, f.normal(), f.represent());

        const Real3 disp_com = direction_com * (l_com / length(direction_com));
        const Real3 disp_ipv = rotate(theta_ipv, f.normal(), dom.ipv());

        // update position
        // XXX to treat polygon surface structure, calculate displacement for
        // each partile without polygon information first
        // and then apply it to each particle with structure

        Particle         p1   = dom.particle_at(0);
        Particle         p2   = dom.particle_at(1);
        const ParticleID pid1 = dom.particle_id_at(0);
        const ParticleID pid2 = dom.particle_id_at(1);

        // ipv is a vector from p1 to p2
        const Real  ratio_p1     = -p1.D() / (p1.D() + p2.D());
        const Real  ratio_p2     =  p2.D() / (p1.D() + p2.D());
        const Real3 disp_ipv_p1  = disp_ipv  * ratio_p1;
        const Real3 disp_ipv_p2  = disp_ipv  * ratio_p2;
        const Real3 disp_ipv0_p1 = dom.ipv() * ratio_p1;
        const Real3 disp_ipv0_p2 = dom.ipv() * ratio_p2;

        Real3 disp_p1 = disp_com - disp_ipv0_p1 + disp_ipv_p1;
        Real3 disp_p2 = disp_com - disp_ipv0_p2 + disp_ipv_p2;
        std::pair<Real3, FaceID> pos_p1(p1.position(), this->get_face_id(pid1));
        std::pair<Real3, FaceID> pos_p2(p2.position(), this->get_face_id(pid2));

        if(0 == ecell4::polygon::travel(this->polygon(), pos_p1, disp_p1, 2))
        {
            SGFRD_TRACE(tracer_.write("escape_com_circular_pair "
                        "p1 moving on face: precision lost"))
        }
        if(0 == ecell4::polygon::travel(this->polygon(), pos_p2, disp_p2, 2))
        {
            SGFRD_TRACE(tracer_.write("escape_com_circular_pair "
                        "p1 moving on face: precision lost"))
        }

        p1.position() = pos_p1.first;
        p2.position() = pos_p2.first;

        this->update_particle(pid1, p1, pos_p1.second);
        this->update_particle(pid2, p2, pos_p2.second);

        boost::array<boost::tuple<ParticleID, Particle, FaceID>, 2> results;
        results[0] = boost::make_tuple(pid1, p1, pos_p1.second);
        results[1] = boost::make_tuple(pid2, p2, pos_p2.second);

        return results;
    }

    boost::container::small_vector<boost::tuple<ParticleID, Particle, FaceID>, 2>
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

    boost::container::small_vector<boost::tuple<ParticleID, Particle, FaceID>, 2>
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

        BOOST_AUTO(const& rules, this->model_->query_reaction_rules(
                    p1.species(), p2.species()));
        assert(false == rules.empty());

        const Real k_tot = this->calc_k_tot(rules);
        boost::optional<ReactionRule const&> optr(boost::none);
        if(rules.size() != 1)
        {
            Real rndr = this->uniform_real() * k_tot;
            BOOST_FOREACH(ReactionRule const& rl, rules)
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

                return boost::container::small_vector<
                    boost::tuple<ParticleID, Particle, FaceID>, 2>(0ul);
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

                const FaceID       sh_fid = sh.structure_id();
                const triangle_type&    f = this->polygon().triangle_at(sh_fid);
                const Real3 direction_com =
                    rotate(theta_com, f.normal(), f.represent());

                Real3 disp_com = direction_com * (l_com / length(direction_com));
                std::pair<Real3, FaceID> pos_com(sh.position(), sh_fid);
                if(0 == ecell4::polygon::travel(
                            this->polygon(), pos_com, disp_com, 2))
                {
                    SGFRD_TRACE(tracer_.write("escape_com_circular_pair "
                                "p1 moving on face: precision lost"))
                }
                const Real3  pos_new = pos_com.first;
                const FaceID fid_new = pos_com.second;

                // make new particle
                const Species species_new =
                    this->model_->apply_species_attributes(rule.products().front());
                const molecule_info_type mol_info =
                    this->world_->get_molecule_info(species_new);
                const Real radius_new = mol_info.radius;
                const Real D_new      = mol_info.D;
                const Particle p_new(species_new, pos_new, radius_new, D_new);

                inside_checker is_inside_of(
                    p_new.position(), p_new.radius(), fid_new, this->polygon());
                if(false == is_inside_of(sh))
                {
                    SGFRD_SCOPE(us, particle_goes_outside, tracer_)
                    const DomainID did = this->get_domain_id(dom);

                    // particle sticks out from the shell.
                    const bool no_overlap =
                        this->burst_and_shrink_overlaps(p_new, fid_new, did);
                    SGFRD_TRACE(tracer_.write("no_overlap = %1%", no_overlap))

                    if(false == no_overlap)
                    {
                        SGFRD_TRACE(tracer_.write(
                                    "reject the reaction because of no space"))
                        const greens_functions::GreensFunction2DRadAbs
                            gf_ipv(dom.D_ipv(), dom.kf(), dom.r0(), dom.sigma(),
                                   dom.R_ipv());
                        const Real theta_ipv = gf_ipv.drawTheta(
                                this->uniform_real(), dom.sigma(), dt);
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
                        if(0 == ecell4::polygon::travel(
                                    this->polygon(), pos_p1, disp_p1, 2))
                        {
                            SGFRD_TRACE(tracer_.write("reaction_ipv_circular_pair "
                                        "p1 moving on face: precision lost"))
                        }
                        if(0 == ecell4::polygon::travel(
                                    this->polygon(), pos_p2, disp_p2, 2))
                        {
                            SGFRD_TRACE(tracer_.write("reaction_ipv_circular_pair "
                                        "p2 moving on face: precision lost"))
                        }
                        p1.position() = pos_p1.first;
                        p2.position() = pos_p2.first;
                        this->update_particle(pid1, p1, pos_p1.second);
                        this->update_particle(pid2, p2, pos_p2.second);

                        boost::container::small_vector<
                            boost::tuple<ParticleID, Particle, FaceID>, 2
                            > results(2);
                        results[0] = boost::make_tuple(pid1, p1, pos_p1.second);
                        results[1] = boost::make_tuple(pid2, p2, pos_p2.second);
                        return results;
                    }
                }
                SGFRD_TRACE(tracer_.write("reaction is accepted."))
                this->update_particle(pid1, p_new, fid_new);

                this->remove_particle(pid2, fid2);
                last_reactions_.push_back(std::make_pair(rule,
                    make_binding_reaction_info(
                        tm, pid1, p1, pid2, p2, pid1, p_new)));

                boost::container::small_vector<
                    boost::tuple<ParticleID, Particle, FaceID>, 2> retval(1ul);
                retval[0] = boost::make_tuple(pid1, p_new, fid_new);
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

        const greens_functions::GreensFunction2DAbsSym
            gf_com(Pair::calc_D_com(p1.D(), p2.D()),
                   Pair::calc_R_com(sh.second.size(), p1, p2));
        const Real t_com_escape = gf_com.drawTime(this->uniform_real());

        const Real k_tot = this->calc_k_tot(this->model_->query_reaction_rules(
                                           p1.species(), p2.species()));
        SGFRD_TRACE(tracer_.write("ipv length = %1%", len_ipv))
        SGFRD_TRACE(tracer_.write("total rate = %1%", k_tot))

        const greens_functions::GreensFunction2DRadAbs
            gf_ipv(Pair::calc_D_ipv(p1.D(), p2.D()),
                   k_tot, len_ipv, p1.radius() + p2.radius(),
                   Pair::calc_R_ipv(sh.second.size(), p1, p2));
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
        SGFRD_TRACE(tracer_.write("particle 1 single reaction time = %1%", t_single_reaction_2))
        SGFRD_TRACE(tracer_.write("single reaction time   = %1%", single_event.first))

        // --------------------------------------------------------------------

        if(gf_event.first < single_event.first)
        {
            SGFRD_TRACE(tracer_.write("gf event occurs first"))
            return Pair(gf_event.second, gf_event.first, this->time(), sh.first,
                        sh.second.size(), std::make_pair(pid1, p1),
                        std::make_pair(pid2, p2), len_ipv, ipv, k_tot);
        }
        else
        {
            SGFRD_TRACE(tracer_.write("single reaction occurs first"))
            return Pair(single_event.second, single_event.first, this->time(),
                        sh.first, sh.second.size(),
                        std::make_pair(pid1, p1), std::make_pair(pid2, p2),
                        len_ipv, ipv, k_tot);
        }
    }


//----------------------------------- multi ------------------------------------

    void fire_multi(Multi& dom, DomainID did)
    {
        SGFRD_SCOPE(us, fire_multi, tracer_);
        SGFRD_TRACE(tracer_.write("fire multi(%1%) for default dt(%2%)", did, dom.dt()));

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
                SGFRD_TRACE(tracer_.write("reaction occurs"))
                std::copy(dom.last_reactions().begin(), dom.last_reactions().end(),
                          std::back_inserter(this->last_reactions_));
                ParticleID pid; Particle p; FaceID fid;
                BOOST_FOREACH(boost::tie(pid, p, fid),
                              this->remove_multi(dom))
                {
                    this->add_event(this->create_closely_fitted_domain(
                        this->create_closely_fitted_shell(pid, p, fid), pid, p));
                }
                SGFRD_TRACE(tracer_.write("multi domain (id = %1%) removed.", did))
                return;
            }
            case Multi::ESCAPE:
            {
                SGFRD_TRACE(tracer_.write("escape occurs"))
                /* burst this domain! */
                ParticleID pid; Particle p; FaceID fid;
                BOOST_FOREACH(boost::tie(pid, p, fid),
                              this->remove_multi(dom))
                {
                    this->add_event(this->create_closely_fitted_domain(
                        this->create_closely_fitted_shell(pid, p, fid), pid, p));
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
        BOOST_FOREACH(boost::tie(pid, p), dom.particles())
        {
            results.push_back(boost::make_tuple(pid, p, this->get_face_id(pid)));
        }

        SGFRD_TRACE(tracer_.write("particles are collected"))

        BOOST_FOREACH(ShellID sid, dom.shell_ids())
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

        BOOST_AUTO(did, get_domain_id(dom));
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

    // XXX: second value of element of result_type is not a mere distance.
    //      - in Multi case, it is just a distance.
    //      - in bursted Single case, it become a distance between centers minus
    //        min_circular_shell of the particle
    std::vector<std::pair<DomainID, Real> >
    burst_and_shrink_non_multis(
            const ParticleID& pid, const Particle& p, const FaceID& fid,
            const std::vector<std::pair<DomainID, Real> >& intruders)
    {
        SGFRD_SCOPE(us, burst_and_shrink_non_multis, tracer_)

        const Real tm(this->time());
        std::vector<std::pair<DomainID, Real> > results;

        DomainID did; Real dist;
        BOOST_FOREACH(boost::tie(did, dist), intruders)
        {
            SGFRD_SCOPE(ns, loop_for_intruders, tracer_)

            BOOST_AUTO(const& ev, get_event(did));
            if(ev->which_domain() == event_type::multi_domain)
            {
                SGFRD_TRACE(tracer_.write("domain %1% is multi", did))
                results.push_back(std::make_pair(did, dist));
                continue;
            }
            SGFRD_TRACE(tracer_.write("domain %1% is not a multi", did))

            DomainID did_; ParticleID pid_; Particle p_; FaceID fid_;
            BOOST_FOREACH(boost::tie(pid_, p_, fid_),
                          burst_event(std::make_pair(did, ev), tm))
            {
                SGFRD_TRACE(tracer_.write(
                    "add closely-fitted domain to bursted particle %1%", pid_))
                did_ = add_event(create_closely_fitted_domain(
                    create_closely_fitted_shell(pid_, p_, fid_), pid_, p_));
                results.push_back(std::make_pair(did_,
                    ecell4::polygon::distance(this->polygon(),
                        std::make_pair(p.position(), fid),
                        std::make_pair(p_.position(), fid_)) -
                    (calc_min_single_circular_shell_radius(p_) * 1.1)));
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
    burst_and_shrink_non_multis(const vertex_id_type& vid,
        const std::vector<std::pair<DomainID, Real> >& intruders)
    {
        SGFRD_SCOPE(us, burst_and_shrink_non_multis_vertex, tracer_)
        const std::pair<Real3, vertex_id_type> vpos = std::make_pair(
                polygon().vertex_at(vid).position, vid);

        const Real tm(this->time());
        std::vector<std::pair<DomainID, Real> > results;

        DomainID did; Real dist;
        BOOST_FOREACH(boost::tie(did, dist), intruders)
        {
            SGFRD_SCOPE(ns, loop_for_intruders, tracer_)

            BOOST_AUTO(const& ev, get_event(did));
            if(ev->which_domain() == event_type::multi_domain)
            {
                SGFRD_TRACE(tracer_.write("domain %1% is multi", did))
                results.push_back(std::make_pair(did, dist));
                continue;
            }
            SGFRD_TRACE(tracer_.write("domain %1% is not a multi", did))

            DomainID did_; ParticleID pid_; Particle p_; FaceID fid_;
            BOOST_FOREACH(boost::tie(pid_, p_, fid_),
                          burst_event(std::make_pair(did, ev), tm))
            {
                SGFRD_TRACE(tracer_.write(
                    "add closely-fitted domain to bursted particle %1%", pid_))

                did_ = add_event(create_closely_fitted_domain(
                    create_closely_fitted_shell(pid_, p_, fid_), pid_, p_));
                results.push_back(std::make_pair(did_,
                    ecell4::polygon::distance(this->polygon(),
                        vpos, std::make_pair(p_.position(), fid_)) -
                    (calc_min_single_circular_shell_radius(p_) * 1.1)));
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

        BOOST_AUTO(const id_of_from, get_domain_id(from));
        remove_event(id_of_from);
        SGFRD_TRACE(tracer_.write("remove domain from(%1%)", id_of_from))

        // reset domain_id
        const domain_id_setter didset(this->get_domain_id(to));
        mut_sh_vis_applier(didset, from);
        SGFRD_TRACE(tracer_.write("set domain ID for shell in from(%1%)", id_of_from))


        // move particle
        ParticleID pid; Particle p; FaceID fid;
        BOOST_FOREACH(boost::tie(pid, p), from.particles())
        {
            const bool addp_result = to.add_particle(pid);
            assert(addp_result);
        }
        // move shell
        BOOST_FOREACH(ShellID sid, from.shell_ids())
        {
            const bool adds_result = to.add_shell(sid);
            SGFRD_TRACE(tracer_.write("Shell(%1%).domain_id = %2%", sid,
                        boost::get<circular_shell_type>(this->get_shell(sid)).domain_id()))
            assert(adds_result);
        }
        to.determine_reaction_length();
        to.determine_delta_t();

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

            escaped_ = false;
            inside_checker is_inside(p.position(), p.radius(), fid, sim.polygon());
            if(applier(is_inside, domain))
            {
                SGFRD_TRACE(sim.tracer_.write("particle is inside of the shell"))
                return true;
            }

            SGFRD_TRACE(sim.tracer_.write("particle escaped"));

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

            SGFRD_TRACE(sim.tracer_.write("particle escaped"))

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

            SGFRD_TRACE(sim.tracer_.write("particle escaped"))

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

        BOOST_AUTO(ev, boost::make_shared<event_type>(dom.begin_time() + dom.dt(), dom));
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
        SGFRDEvent::domain_type& dom = ev.second->domain();

        switch(ev.second->which_domain())
        {
            case event_type::single_domain:
            {
                return this->fire_single(boost::get<Single>(dom), ev.first);
            }
            case event_type::pair_domain:
            {
                return this->fire_pair(boost::get<Pair>(dom), ev.first);
            }
            case event_type::multi_domain:
            {
                return this->fire_multi(boost::get<Multi>(dom), ev.first);
            }
            case event_type::birth_domain:
            {
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

    ShellID create_closely_fitted_shell(
            const ParticleID& pid, const Particle& p, const FaceID fid)
    {
        SGFRD_SCOPE(us, create_closely_fitted_shell, tracer_);
        SGFRD_TRACE(tracer_.write("add close shell for %1% @ face %2%", pid, fid))
        SGFRD_TRACE(tracer_.write("species %1% has radius %2%",
                    p.species_serial(), p.radius()))

        const ShellID sid(shell_id_gen());
        circular_shell_type sh(circle_type(p.radius(), p.position(),
                               this->polygon().triangle_at(fid).normal()), fid);
        SGFRD_TRACE(tracer_.write("shell has size == %1%", sh.size()))

        shell_container_.add_shell(sid, sh, fid);
        SGFRD_TRACE(tracer_.write("new shell id is %1%", sid))
        return sid;
    }
    Single create_closely_fitted_domain(
            const ShellID& sid, const ParticleID& pid, const Particle& p)
    {
        SGFRD_SCOPE(us, create_closely_fitted_domain, tracer_);
        SGFRD_TRACE(tracer_.write("for particle %1%, shell %2%", pid, sid))
        return Single(Single::ESCAPE, 0., this->time(), sid, std::make_pair(pid, p));
    }

    std::pair<ShellID, circle_type>
    create_minimum_single_shell(
            const ParticleID& pid, const Particle& p, const FaceID fid)
    {
        return create_single_circular_shell(std::make_pair(p.position(), fid),
            calc_min_single_circular_shell_radius(p));
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


    std::vector<std::pair<vertex_id_type, Real> >
    get_intrusive_vertices(const std::pair<Real3, FaceID>& pos,
                           const Real radius) const
    {
        return polygon().list_vertices_within_radius(pos, radius);
    }

    std::vector<std::pair<DomainID, Real> >
    get_intrusive_domains(const std::pair<Real3, FaceID>& pos,
                          const Real radius) const
    {
        SGFRD_SCOPE(us, get_intrusive_domains_position, tracer_);

        const std::vector<std::pair<std::pair<ShellID, shell_type>, Real>
            > shells(shell_container_.list_shells_within_radius(pos, radius));

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
            boost::tie(shell_id_pair, dist) = *iter;

            SGFRD_TRACE(tracer_.write("shell %1% = {%2%} is at %3% distant",
                        shell_id_pair.first, shell_id_pair.second, dist));

            const DomainID did = boost::apply_visitor(
                    domain_id_getter(), shell_id_pair.second);

            SGFRD_TRACE(tracer_.write("shell %1% is related to domain %2%",
                        shell_id_pair.first, did));

            // filter by domainID to be unique
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

    std::vector<std::pair<DomainID, Real> >
    get_intrusive_domains(const vertex_id_type& vid, const Real radius) const
    {
        SGFRD_SCOPE(us, get_intrusive_domains_vid, tracer_);

        const std::pair<Real3, vertex_id_type> vpos = std::make_pair(
                polygon().vertex_at(vid).position, vid);
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
            boost::tie(shell_id_pair, dist) = *iter;
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
        Real lensq = std::numeric_limits<Real>::max();
        const boost::array<ecell4::Segment, 6>& barrier =
            polygon().face_at(pos.second).barrier;

        for(std::size_t i=0; i<6; ++i)
        {
            const Real dist2 = ecell4::sgfrd::distance_sq(pos.first, barrier[i]);
            if(dist2 < lensq)
            {
                lensq = dist2;
            }
        }
        return std::sqrt(lensq);
    }
    Real get_max_cone_size(const vertex_id_type& vid) const
    {
        return polygon().vertex_at(vid).max_conical_shell_size * 0.5;
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
        BOOST_FOREACH(ReactionRule const& rule, rules)
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
        BOOST_FOREACH(ReactionRule const& rule, rules)
        {
            a += rule.k();
            if(a > thrs) return rule;
        }
        throw std::logic_error("reaction cannot detemined");
    }

    // for BD part...
    Real3 random_circular_uniform(const Real r, const FaceID& fid)
    {
        const Real theta = this->rng_.uniform(0., 2 * M_PI);
        const Real3 rnd(r * std::cos(theta), r * std::sin(theta), 0.);
        const Real3& normal = this->polygon().triangle_at(fid).normal();
        const Real tilt  = angle(Real3(0, 0, 1), normal);

             if(std::abs(tilt - M_PI) < 1e-8){return rnd;}
        else if(std::abs(tilt + M_PI) < 1e-8){return rnd * (-1.0);}

        const Real3 axis = cross_product(Real3(0., 0., 1.), normal);
        return rotate(tilt, axis * (1. / length(axis)), rnd);
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

};


} // sgfrd
} // ecell4
#endif // ECELL4_SGFRD_SIMULATOR
