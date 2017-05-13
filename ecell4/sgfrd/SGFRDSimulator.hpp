#ifndef ECELL4_SGFRD_SIMULATOR
#define ECELL4_SGFRD_SIMULATOR

#include <greens_functions/GreensFunction2DAbsSym.hpp>
#include <greens_functions/GreensFunction2DRefWedgeAbs.hpp>
#include <ecell4/core/SimulatorBase.hpp>
#include <ecell4/core/SerialIDGenerator.hpp>
#include <ecell4/core/geometry.hpp>
#include "ShellContainer.hpp"
#include "ShellVisitorApplier.hpp"
#include "ShellVisitors.hpp"
#include "SGFRDEvent.hpp"
#include "SGFRDWorld.hpp"
#include <boost/tuple/tuple.hpp>
#include <iostream>

#ifndef SGFRD_NDEBUG
#define DUMP_MESSAGE( str ) std::cerr << str << std::endl
#else
#define DUMP_MESSAGE( str )
#endif //NDEBUG

namespace ecell4
{
namespace sgfrd
{

template<typename T_polygon_traits>
class SGFRDSimulator :
    public ecell4::SimulatorBase<ecell4::Model, SGFRDWorld<T_polygon_traits> >
{
  public:
    // polygon
    typedef T_polygon_traits    polygon_traits_type;
    typedef ecell4::Polygon<polygon_traits_type>     polygon_type;
    typedef typename polygon_type::triangle_type     triangle_type;
    typedef typename polygon_type::face_id_type      face_id_type;
    typedef typename polygon_type::edge_id_type      edge_id_type;
    typedef typename polygon_type::vertex_id_type    vertex_id_type;
    typedef typename polygon_type::face_descripter   face_descripter;
    typedef typename polygon_type::edge_descripter   edge_descripter;
    typedef typename polygon_type::vertex_descripter vertex_descripter;
    typedef typename polygon_type::local_index_type  local_index_type;
    typedef typename polygon_type::barycentric_type  barycentric_type;

    // Event
    typedef SGFRDEvent          event_type;
    typedef SGFRDEventScheduler scheduler_type;
    typedef typename SGFRDEvent::domain_type domain_type;
    typedef typename SGFRDEventScheduler::value_type event_id_pair_type;

    // Simulator
    typedef ecell4::SimulatorBase<ecell4::Model, SGFRDWorld<polygon_traits_type>
            > base_type;
    typedef typename base_type::world_type world_type;
    typedef typename base_type::model_type model_type;

    // ShellContainer
    typedef ecell4::SerialIDGenerator<ShellID> shell_id_generator_type;
    typedef ShellContainer<polygon_traits_type> shell_container_type;
    typedef typename shell_container_type::shell_type           shell_type;
    typedef typename shell_container_type::shell_id_pair_type   shell_id_pair_type;
    typedef typename shell_container_type::circle_type          circle_type;
    typedef typename shell_container_type::conical_surface_type conical_surface_type;
    typedef typename shell_container_type::circular_shell_type  circular_shell_type;
    typedef typename shell_container_type::conical_surface_shell_type
            conical_surface_shell_type;
    typedef mutable_shell_visitor_applier<polygon_traits_type>
            mutable_shell_visitor_applier_type;

  public:

    SGFRDSimulator(const boost::shared_ptr<world_type>& world,
                   const boost::shared_ptr<model_type>& model,
                   Real bd_dt_factor = 1e-5)
        : base_type(model, world), dt_(0), bd_dt_factor_(bd_dt_factor),
          rng_(*(world->rng())), shell_container_(world->polygon()),
          mut_sh_vis_applier(shell_container_)
    {}

    SGFRDSimulator(boost::shared_ptr<world_type> world, Real bd_dt_factor = 1e-5)
        : base_type(world), dt_(0), bd_dt_factor_(bd_dt_factor),
          rng_(*(world->rng())), shell_container_(world->polygon()),
          mut_sh_vis_applier(shell_container_)
    {}

    ~SGFRDSimulator(){}

    void step();
    bool step(const Real& upto);

    void initialize();
    void finalize();

    Real dt() const {return dt_;}

//     bool check_reaction() const {return last_reactions_.size() > 0;}
//     std::vector<std::pair<ReactionRule, reaction_info_type> >
//     last_reactions() const {return last_reactions_;}

  private: // wrappers

    Real uniform_real(){return this->rng_.random();}

    world_type   const& world()   const {return *(this->world_);}
    polygon_type const& polygon() const {return this->world_->polygon();}

    bool update_particle(const ParticleID& pid, const Particle& p,
                         const face_id_type& fid)
    {return this->world_->update_particle(pid, p, fid);}


    shell_type&       get_shell(ShellID const& id)
    {return shell_container_.get_shell(id);}
    shell_type const& get_shell(ShellID const& id) const
    {return shell_container_.get_shell(id);}
    void remove_shell(ShellID const& id)
    {return shell_container_.remove_shell(id);}

    face_id_type      get_face_id(const ParticleID& pid) const
    {return this->world_->get_face_id(pid);}
    std::pair<ParticleID, Particle> get_particle(const ParticleID& pid) const
    {return this->world_->get_particle(pid);}

    // consider scheduler.time
    Real time() const {return this->world_->t();}
    void set_time(const Real t) {return this->world_->set_t(t);}

  private:

    //! make event from domain and push it into scheduler
    template<typename domainT>
    void add_event(const domainT& dom);
    //! fire the event
    void fire_event(const event_id_pair_type& event);
    struct domain_firer;

    //! burst the event and store the bursted domain to result
    void burst_event(const event_type& event,
                     std::vector<domain_type>& result);
    struct domain_burster;
    struct single_shell_burster;
    struct pair_shell_burster;

    void burst_non_multis(std::vector<domain_type>& domains);

    //! single event. after fired, it creates next event.
    struct single_shell_escapement;
    struct single_shell_reactor;
    void escape  (const Single& domain);
    void reaction(const Single& domain);

    //! make shell that has same size as particle radius. add the shell to scon.
    ShellID create_minimum_shell(
            const ParticleID& pid, const Particle& p, const face_id_type fid);
    //! make single that has same shell size as particle radius.
    Single create_minimum_domain(
            const ShellID& sid, const ParticleID& pid, const Particle& p);

    //! create single circular shell from geometric & domain information
    //  this adds the shell to shell container.
    std::pair<ShellID, circle_type>
    create_single_circular_shell(
            const std::pair<Real3, face_id_type>& pos, const Real size);

    //! create single conical surf shell from geometric & domain information
    //  this adds the shell to shell container.
    std::pair<ShellID, conical_surface_type>
    create_single_conical_surface_shell(
            const vertex_id_type& pos, const Real size);

    //! create domain and determine EventKind using GF
    Single create_single(const std::pair<ShellID, circle_type>& sh,
            const ParticleID& pid, const Particle& p);
    Single create_single(const std::pair<ShellID, conical_surface_type>& sh,
            const ParticleID& pid, const Particle& p);

    //! make domain and call add_event
    void create_event(
            const ParticleID& pid, const Particle& p, const face_id_type fid);


    std::vector<std::pair<vertex_id_type, Real> >
    get_intrusive_vertices(const std::pair<Real3, face_id_type>& pos,
                           const Real radius) const;

    std::vector<std::pair<DomainID, Real> >
    get_intrusive_domains(const std::pair<Real3, face_id_type>& pos,
                          const Real radius) const;

    std::vector<std::pair<DomainID, Real> >
    get_intrusive_domains(const vertex_id_type& vid, const Real radius) const;

    //! calculate max circle size allowed by geometric restraints
    Real get_max_circle_size(const std::pair<Real3, face_id_type>& pos) const;
    //! calculate max cone size allowed by geometric restraints
    Real get_max_cone_size(const vertex_id_type& vid) const;

  private:

    static const Real single_circular_shell_factor;
    static const Real single_conical_surface_shell_factor;

  private:

    // from SimulatorBase
    // boost::shared_ptr<model_type> model_;
    // boost::shared_ptr<world_type> world_;
    // Integer num_steps_;
    Real dt_;
    Real bd_dt_factor_;
    ecell4::RandomNumberGenerator&     rng_;
    scheduler_type                     scheduler_;
    shell_id_generator_type            shell_id_gen;
    shell_container_type               shell_container_;
    mutable_shell_visitor_applier_type mut_sh_vis_applier;
//     std::vector<std::pair<ReactionRule, reaction_info_type> > last_reactions_;
};

template<typename T>
const Real SGFRDSimulator<T>::single_circular_shell_factor = 1.5;
template<typename T>
const Real SGFRDSimulator<T>::single_conical_surface_shell_factor = 1.5;

template<typename T>
void SGFRDSimulator<T>::step()
{
    this->set_time(this->scheduler_.next_time());
    DUMP_MESSAGE("start firing event");
    fire_event(this->scheduler_.pop());
    DUMP_MESSAGE("now " << shell_container_.num_shells() << " shells exist.");
    return;
}

template<typename T>
bool SGFRDSimulator<T>::step(const Real& upto)
{
    this->step();
    return this->time() < upto;
}

template<typename T>
void SGFRDSimulator<T>::initialize()
{
    std::vector<std::pair<ParticleID, Particle> > const& ps =
        this->world_->list_particles();
    DUMP_MESSAGE("number of particles: " << ps.size());
    for(std::vector<std::pair<ParticleID, Particle> >::const_iterator
        iter = ps.begin(); iter != ps.end(); ++iter)
    {
        DUMP_MESSAGE("making minimum shells...");
        add_event(create_minimum_domain(create_minimum_shell(
                      iter->first, iter->second, this->get_face_id(iter->first)),
                  iter->first, iter->second));
    }
    return ;
}

template<typename T>
void SGFRDSimulator<T>::finalize()
{
    std::vector<domain_type> tmp;
    while(scheduler_.size() != 0)
    {
        burst_event(*(this->scheduler_.pop().second), tmp);
    }
    return ;
}

template<typename T>
template<typename domainT>
void SGFRDSimulator<T>::add_event(const domainT& dom)
{
    const DomainID did = scheduler_.add(
            boost::make_shared<event_type>(dom.begin_time() + dom.dt(), dom));
    domain_id_setter didset(did);
    mut_sh_vis_applier(didset, dom);
    return ;
}

template<typename T>
struct SGFRDSimulator<T>::domain_firer : boost::static_visitor<void>
{
    domain_firer(SGFRDSimulator<T>& s) : sim(s){}

    void operator()(const Single& dom)
    {
        switch(dom.eventkind())
        {
            case Single::ESCAPE:
            {
                DUMP_MESSAGE("single escape");
                sim.escape(dom);
                return;
            }
            case Single::REACTION:
            {
                DUMP_MESSAGE("single reaction");
                sim.reaction(dom);
                return;
            }
            case Single::UNKNOWN:
                throw std::logic_error("when firing Single: event unspecified");
            default:
                throw std::logic_error("when firing Single: invalid enum value");
        }
    }

    void operator()(const Pair& dom)
    {//TODO
        std::cerr << "[WARNING] Pair domain firer has not been implemented yet"
                  << std::endl;
        return;
    }

    void operator()(const Multi& dom)
    {//TODO
        std::cerr << "[WARNING] Multi domain firer has not been implemented yet"
                  << std::endl;
        return;
    }

  private:

    SGFRDSimulator<T>& sim;
};

template<typename T>
void SGFRDSimulator<T>::fire_event(const event_id_pair_type& ev)
{
    DUMP_MESSAGE("applying firer");
    domain_firer firer(*this);
    boost::apply_visitor(firer, ev.second->domain());
    DUMP_MESSAGE("end applying firer");
    return ;
}

template<typename T>
struct SGFRDSimulator<T>::single_shell_burster : boost::static_visitor<void>
{
    single_shell_burster(
            SGFRDSimulator<T>& s, const domain_type& d, std::vector<domain_type>& r)
        : sim(s), dom(d), result(r)
    {}

    void operator()(const circular_shell_type& sh)
    {//TODO
        std::cerr << "bursting circualr shell" << std::endl;
        return;
    }

    void operator()(const conical_surface_shell_type& sh)
    {//TODO
        std::cerr << "bursting conical surface shell" << std::endl;
        return;
    }

  private:

    SGFRDSimulator<T>& sim;
    const domain_type& dom;
    std::vector<domain_type>& result;
};

template<typename T>
struct SGFRDSimulator<T>::pair_shell_burster : boost::static_visitor<void>
{
    pair_shell_burster(SGFRDSimulator<T>& s, const domain_type& d,
                       std::vector<domain_type>& r)
        : sim(s), dom(d), result(r)
    {}

    void operator()(const circular_shell_type& sh)
    {//TODO
        std::cerr << "bursting pair circualr shell" << std::endl;
        return;
    }

    void operator()(const conical_surface_shell_type& sh)
    {
        throw std::logic_error("bursting pair conical surface shell");
    }

  private:

    SGFRDSimulator<T>& sim;
    const domain_type& dom;
    std::vector<domain_type>& result;
};


template<typename T>
struct SGFRDSimulator<T>::domain_burster : boost::static_visitor<void>
{
    domain_burster(SGFRDSimulator<T>& s, std::vector<domain_type>& r)
        : sim(s), result(r)
    {}

    void operator()(const Single& dom)
    {
        single_shell_burster burster(sim, dom, result);
        boost::apply_visitor(burster, sim.get_shell(dom.shell_id()));
        sim.remove_shell(dom.shell_id());
        return;
    }

    void operator()(const Pair& dom)
    {
        pair_shell_burster burster(sim, dom, result);
        boost::apply_visitor(burster, sim.get_shell(dom.shell_id()));
        sim.remove_shell(dom.shell_id());
        return;
    }

    void operator()(const Multi& dom)
    {
        throw std::logic_error("bursting Multi Shell");
    }

  private:

    SGFRDSimulator<T>& sim;
    std::vector<domain_type>& result;
};

template<typename T>
void SGFRDSimulator<T>::burst_event(
        const event_type& ev, std::vector<domain_type>& res)
{
    domain_burster burster(*this, res);
    boost::apply_visitor(burster, ev.domain());
    return ;
}

template<typename T>
struct SGFRDSimulator<T>::single_shell_escapement : boost::static_visitor<void>
{
    single_shell_escapement(SGFRDSimulator& s, const Single& d)
        : sim(s), dom(d)
    {}

    void operator()(const circular_shell_type& sh)
    {
        DUMP_MESSAGE("single shell escapement circular shell");
        if(sh.size() == dom.particle().radius())
        {
            DUMP_MESSAGE("minimum shell. didnot move.");
            results = boost::make_tuple(dom.particle_id(), dom.particle(),
                                        sim.get_face_id(dom.particle_id()));
            return;
        }

        Particle   p   = dom.particle();
        ParticleID pid = dom.particle_id();

        const Real r   = sh.size() - p.radius();
        const Real theta = sim.uniform_real() * 2.0 * 3.141592653589793;
        DUMP_MESSAGE("r = " << r << ", theta = " << theta);
        const face_id_type   fid  = sim.get_face_id(pid);
        const triangle_type& face = sim.polygon().triangle_at(fid);
        const Real3 direction = rotate(theta, face.normal(), face.represent());
        DUMP_MESSAGE("direction = " << direction << ", length = " << length(direction));

        std::pair<std::pair<Real3, face_id_type>, Real3> state =
            std::make_pair(/*position = */std::make_pair(p.position(), fid),
                       /*displacement = */direction * r / length(direction));

        DUMP_MESSAGE("pos  = " << state.first.first << ", fid = " << state.first.second);
        unsigned int continue_count = 2;
        while(continue_count > 0)
        {
            state = sim.polygon().move_next_face(state.first, state.second);
            const Real3& disp = state.second;
            if(disp[0] == 0. && disp[1] == 0. && disp[2] == 0.) break;
            --continue_count;
            DUMP_MESSAGE("pos  = " << state.first.first << ", fid = " << state.first.second);
            DUMP_MESSAGE("disp = " << disp << ", length = " << length(disp));
        }
        if(continue_count == 0)
            std::cerr << "[WARNING] moving on face: precision lost" << std::endl;

        DUMP_MESSAGE("pos  = " << state.first.first << ", fid = " << state.first.second);
        DUMP_MESSAGE("disp = " << state.second << ", length = " << state.second);

        DUMP_MESSAGE("escaped.");

        p.position() = state.first.first;
        sim.update_particle(pid, p, state.first.second);
        results = boost::make_tuple(pid, p, state.first.second);
        return ;
    }

    void operator()(const conical_surface_shell_type& sh)
    {
        Particle           p   = dom.particle();
        const ParticleID   pid = dom.particle_id();
        const face_id_type fid = sim.get_face_id(pid);
        DUMP_MESSAGE("escape-conical: pos  = " << p.position() << ", fid = " << fid);

        const Real r     = sh.size() - p.radius();
        greens_functions::GreensFunction2DRefWedgeAbs
            gf(p.D(), length(p.position() - sh.position()),
               r,     sh.shape().apex_angle());
        const Real theta = gf.drawTheta(sim.uniform_real(), r, dom.dt());

        DUMP_MESSAGE("escape-conical: r = " << r << ", theta = " << theta);

        const std::pair<Real3, face_id_type> state =
            sim.polygon().rotate_around_vertex(std::make_pair(p.position(), fid),
                                               sh.structure_id(), r, theta);

        DUMP_MESSAGE("escaped : pos = " << state.first << ", fid = " << state.second);

        p.position() = state.first;
        sim.update_particle(pid, p, state.second);
        results = boost::make_tuple(pid, p, state.second);
        return;
    }

    // single shell escapement result can be only one particle
    boost::tuple<ParticleID, Particle, face_id_type> results;

  private:

    SGFRDSimulator<T>& sim;
    Single const& dom;

};

template<typename T>
void SGFRDSimulator<T>::escape(const Single& domain)
{
    DUMP_MESSAGE("single escapement begin");
    single_shell_escapement escapement(*this, domain);
    boost::apply_visitor(escapement, get_shell(domain.shell_id()));
    DUMP_MESSAGE("single escapement applied");
    remove_shell(domain.shell_id());
    this->create_event(boost::get<0>(escapement.results),
                       boost::get<1>(escapement.results),
                       boost::get<2>(escapement.results));
    DUMP_MESSAGE("shell removed");
    return;
}

template<typename T>
struct SGFRDSimulator<T>::single_shell_reactor : boost::static_visitor<void>
{
    void operator()(const circular_shell_type& sh)
    {//TODO
        std::cerr << "react in single circualr shell" << std::endl;
        return;
    }

    void operator()(const conical_surface_shell_type& sh)
    {//TODO
        std::cerr << "react in single conical surface shell" << std::endl;
        return;
    }
};

template<typename T>
void SGFRDSimulator<T>::reaction(const Single& domain)
{
    single_shell_reactor reactor;
    boost::apply_visitor(reactor, get_shell(domain.shell_id()));
    remove_shell(domain.shell_id());
    return;
}

template<typename T>
ShellID SGFRDSimulator<T>::create_minimum_shell(
        const ParticleID& pid, const Particle& p, const face_id_type fid)
{
    const ShellID sid(shell_id_gen());
    circular_shell_type sh(
            circle_type(p.radius(), p.position(),
                        this->polygon().triangle_at(fid).normal()),
            fid);

    shell_container_.add_shell(sid, sh, fid);
    return sid;
}

template<typename T>
Single SGFRDSimulator<T>::create_minimum_domain(const ShellID& sid,
        const ParticleID& pid, const Particle& p)
{
    return Single(Single::ESCAPE, 0., this->time(), sid, std::make_pair(pid, p));
}

template<typename T>
void SGFRDSimulator<T>::create_event(
            const ParticleID& pid, const Particle& p, const face_id_type fid)
{
    DUMP_MESSAGE("create event");
    const std::pair<Real3, face_id_type> pos = std::make_pair(p.position(), fid);

    const Real min_circle_size = p.radius() * single_circular_shell_factor;
          Real max_circle_size = get_max_circle_size(pos);

    DUMP_MESSAGE("min circle size = " << min_circle_size);
    DUMP_MESSAGE("max circle size = " << max_circle_size);

    const Real vertices_search_range = (max_circle_size > min_circle_size) ? max_circle_size :
        std::numeric_limits<Real>::infinity();
    DUMP_MESSAGE("vertices range = " << vertices_search_range);

    const std::vector<std::pair<vertex_id_type, Real> > intrusive_vertices(
            get_intrusive_vertices(pos, vertices_search_range));

    const bool draw_circular = max_circle_size < min_circle_size ? false :
                               intrusive_vertices.empty() ? true :
                              (intrusive_vertices.front().second > min_circle_size);
    if(draw_circular)
    {
        DUMP_MESSAGE("drawing circular shell");
        if(!intrusive_vertices.empty())
        {
            DUMP_MESSAGE("vertices found. " << intrusive_vertices.front().second);
            max_circle_size = intrusive_vertices.front().second;
        }

        DUMP_MESSAGE("creating single circle that size is " << max_circle_size);

        //XXX! consider shell-removing timing !
        const std::vector<std::pair<DomainID, Real> > intrusive_domains(
                get_intrusive_domains(pos, max_circle_size));

        DUMP_MESSAGE("intrusive domains: " << intrusive_domains.size());

        if(intrusive_domains.empty())
        {
            return add_event(create_single(create_single_circular_shell(
                             pos, max_circle_size), pid, p));
        }

        if(intrusive_domains.front().second > min_circle_size)
        {
            return add_event(create_single(create_single_circular_shell(
                             pos, intrusive_domains.front().second), pid, p));
        }

        // TODO burst!
        std::cerr << "[WARNING] ignoring intrusive domains." << std::endl;
        return add_event(create_single(create_single_circular_shell(
                         pos, max_circle_size), pid, p));
    }
    else // conical surface shell
    {
        DUMP_MESSAGE("drawing conical shell");
        const vertex_id_type& vid = intrusive_vertices.front().first;
        DUMP_MESSAGE("vertex id = " << vid << ", distance = " << intrusive_vertices.front().second);
        const Real min_cone_size = p.radius() * single_conical_surface_shell_factor;
        const Real max_cone_size = get_max_cone_size(vid);
        DUMP_MESSAGE("min cone size = " << min_cone_size);
        DUMP_MESSAGE("max cone size = " << max_cone_size);
        const std::vector<std::pair<DomainID, Real> > intrusive_domains(
                get_intrusive_domains(vid, max_cone_size));
        DUMP_MESSAGE("intrusive domains = " << intrusive_domains.size());

        if(intrusive_domains.empty())
        {
            return add_event(create_single(create_single_conical_surface_shell(
                                           vid, max_cone_size), pid, p));
        }

        if(intrusive_domains.front().second > min_cone_size)
        {
            DUMP_MESSAGE("avoid domains overlapping");
            return add_event(create_single(create_single_conical_surface_shell(
                             vid, intrusive_domains.front().second), pid, p));
        }

        // TODO burst and form pair or multi if needed
        std::cerr << "[WARNING] intrusive domains exist." << std::endl;
        return add_event(create_single(create_single_conical_surface_shell(
                                       vid, max_cone_size), pid, p));
    }
}

template<typename T>
std::pair<ShellID, typename SGFRDSimulator<T>::circle_type>
SGFRDSimulator<T>::create_single_circular_shell(
            const std::pair<Real3, face_id_type>& pos, const Real size)
{
    DUMP_MESSAGE("create single circular shell");
    const ShellID id(shell_id_gen());
    const circle_type shape(size, pos.first,
                            polygon().triangle_at(pos.second).normal());
    shell_container_.add_shell(id, circular_shell_type(shape, pos.second),
                               pos.second);
    return std::make_pair(id, shape);
}

template<typename T>
std::pair<ShellID, typename SGFRDSimulator<T>::conical_surface_type>
SGFRDSimulator<T>::create_single_conical_surface_shell(
            const vertex_id_type& vid, const Real size)
{
    DUMP_MESSAGE("create single conical surface shell");
    const ShellID id(shell_id_gen());
    const conical_surface_type shape(polygon().vertex_at(vid).position,
                                     polygon().apex_angle(vid), size);

    shell_container_.add_shell(id, conical_surface_shell_type(shape, vid),
                               vid);
    return std::make_pair(id, shape);
}


template<typename T>
Single SGFRDSimulator<T>::create_single(
        const std::pair<ShellID, circle_type>& sh,
        const ParticleID& pid, const Particle& p)
{
    //TODO consider single-reaction
    DUMP_MESSAGE("create single domain having circular shell");

    const greens_functions::GreensFunction2DAbsSym
        gf(/* D = */ p.D(),
           /* a = */ sh.second.size() - p.radius());
    const Real dt = gf.drawTime(uniform_real());
    DUMP_MESSAGE("delta t calculated: " << dt);

    return Single(Single::ESCAPE, dt, this->time(), sh.first,
                  std::make_pair(pid, p));
}

template<typename T>
Single SGFRDSimulator<T>::create_single(
        const std::pair<ShellID, conical_surface_type>& sh,
        const ParticleID& pid, const Particle& p)
{
    //TODO consider single-reaction
    DUMP_MESSAGE("create single domain having conical shell");
    DUMP_MESSAGE("D   = " << p.D());
    DUMP_MESSAGE("r0  = " << length(p.position() - sh.second.apex()));
    DUMP_MESSAGE("a   = " << sh.second.size() - p.radius());
    DUMP_MESSAGE("phi = " << sh.second.apex_angle());

    const greens_functions::GreensFunction2DRefWedgeAbs
        gf(/* D   = */ p.D(),
           /* r0  = */ length(p.position() - sh.second.apex()),
           /* a   = */ sh.second.size() - p.radius(),
           /* phi = */ sh.second.apex_angle());
    const Real dt = gf.drawTime(uniform_real());

    return Single(Single::ESCAPE, dt, this->time(), sh.first,
                  std::make_pair(pid, p));
}

template<typename T>
std::vector<std::pair<typename SGFRDSimulator<T>::vertex_id_type, Real> >
SGFRDSimulator<T>::get_intrusive_vertices(
        const std::pair<Real3, face_id_type>& pos, const Real radius) const
{
    return polygon().list_vertices_within_radius(pos, radius);
}

template<typename T>
std::vector<std::pair<DomainID, Real> >
SGFRDSimulator<T>::get_intrusive_domains(
        const std::pair<Real3, face_id_type>& pos, const Real radius) const
{
    const std::vector<std::pair<std::pair<ShellID, shell_type>, Real>
        > shells(shell_container_.list_shells_within_radius(pos, radius));

    std::vector<std::pair<DomainID, Real> > domains(shells.size());

    typename std::vector<std::pair<DomainID, Real> >::iterator
        dest = domains.begin();
    for(typename std::vector<std::pair<std::pair<ShellID, shell_type>, Real>
        >::const_iterator iter(shells.begin()), end(shells.end()); iter != end; ++iter)
    {
        *dest = std::make_pair(
                boost::apply_visitor(domain_id_getter(), iter->first.second),
                iter->second);
        ++dest;
    }
    return domains;
}

template<typename T>
std::vector<std::pair<DomainID, Real> >
SGFRDSimulator<T>::get_intrusive_domains(
        const vertex_id_type& vid, const Real radius) const
{
    const std::pair<Real3, vertex_id_type> vpos = std::make_pair(
            polygon().vertex_at(vid).position, vid);
    const std::vector<std::pair<std::pair<ShellID, shell_type>, Real>
        > shells(shell_container_.list_shells_within_radius(vpos, radius));

    std::vector<std::pair<DomainID, Real> > domains(shells.size());

    typename std::vector<std::pair<DomainID, Real> >::iterator
        dest = domains.begin();
    for(typename std::vector<std::pair<std::pair<ShellID, shell_type>, Real>
        >::const_iterator iter(shells.begin()), end(shells.end()); iter != end; ++iter)
    {
        *dest = std::make_pair(
                boost::apply_visitor(domain_id_getter(), iter->first.second),
                iter->second);
        ++dest;
    }
    return domains;
}

template<typename T>
Real SGFRDSimulator<T>::get_max_circle_size(
        const std::pair<Real3, face_id_type>& pos) const
{
    Real lensq = std::numeric_limits<Real>::max();
    const boost::array<std::pair<Real3, Real3>, 6>& barrier =
        polygon().face_at(pos.second).barrier;

    for(std::size_t i=0; i<6; ++i)
    {
        //XXX distance to the segment!
        const Real3 a = pos.first         - barrier[i].first;
        const Real3 b = barrier[i].second - barrier[i].first;
        const Real dot = dot_product(a, b);
        const Real dist2 = length_sq(a) - dot * dot / length_sq(b);

        if(lensq > dist2)
        {
            DUMP_MESSAGE("max_circle_size updated: " << i << ", " << std::sqrt(dist2));
            DUMP_MESSAGE("barrier: " << barrier[i].first[0]  << " " << barrier[i].first[1]  << " " << barrier[i].first[2]);
            DUMP_MESSAGE("barrier: " << barrier[i].second[0] << " " << barrier[i].second[1] << " " << barrier[i].second[2]);
            lensq = dist2;
        }
    }
    return std::sqrt(lensq);
}

template<typename T>
Real SGFRDSimulator<T>::get_max_cone_size(const vertex_id_type& vid) const
{
    return polygon().vertex_at(vid).max_conical_shell_size * 0.5;
}

} // sgfrd
} // ecell4
#endif // ECELL4_SGFRD_SIMULATOR
