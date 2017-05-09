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
#include <iostream>

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
    typedef DomainID            domain_id_type;
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
          shell_container_(world->polygon()), mut_sh_vis_applier(shell_container_)
    {
        initialize();
    }

    SGFRDSimulator(boost::shared_ptr<world_type> world, Real bd_dt_factor = 1e-5)
        : base_type(world), dt_(0), bd_dt_factor_(bd_dt_factor),
          shell_container_(world->polygon()), mut_sh_vis_applier(shell_container_)
    {
        initialize();
    }

    ~SGFRDSimulator(){}

    void step();
    bool step(const Real& upto);

    void initialize();
    void finalize();

//     bool check_reaction() const {return last_reactions_.size() > 0;}
//     std::vector<std::pair<ReactionRule, reaction_info_type> >
//     last_reactions() const {return last_reactions_;}

    boost::shared_ptr<RandomNumberGenerator> rng() {return this->world_->rng();}

  private: // wrappers

    world_type   const& world()   const {return *(this->world_);}
    polygon_type const& polygon() const {return this->world_->polygon();}

    bool update_particle(const ParticleID& pid, const Particle& p,
                         const face_id_type& fid)
    {return this->world_->update_particle(pid, p);}

    Real uniform_real() {this->world_->rng()->random();}

    shell_type&       get_shell(ShellID const& id)
    {return shell_container_.get_shell(id);}
    shell_type const& get_shell(ShellID const& id) const
    {return shell_container_.get_shell(id);}
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

    //! make shell that has same size as particle radius
    ShellID create_minimum_shell(
            const ParticleID& pid, const Particle& p, const face_id_type fid);
    //! make single that has same shell size as particle radius
    Single create_minimum_domain(
            const ShellID& sid, const ParticleID& pid, const Particle& p);

    //! create single circular shell from geometric & domain information
    std::pair<ShellID, circle_type>
    create_single_circular_shell(
            const std::pair<Real3, face_id_type>& pos, const Real size);

    //! create single conical surf shell from geometric & domain information
    std::pair<ShellID, conical_surface_type>
    create_single_conical_surface_shell(
            const std::pair<Real3, vertex_id_type>& pos, const Real size);

    //! create domain and determine EventKind using GF
    Single create_single(const std::pair<ShellID, circle_type>& sh,
            const ParticleID& pid, const Particle& p);
    Single create_single(const std::pair<ShellID, conical_surface_type>& sh,
            const ParticleID& pid, const Particle& p);

    //! make domain and call add_event
    void create_event(
            const ParticleID& pid, const Particle& p, const face_id_type fid);


    std::pair<std::vector<vertex_id_type>, std::pair<vertex_id_type, Real> >
    get_intrusive_vertices(const std::pair<Real3, face_id_type>& pos,
                           const Real radius) const;

    std::pair<std::vector<domain_id_type>, std::pair<domain_id_type, Real> >
    get_intrusive_domains(const std::pair<Real3, face_id_type>& pos,
                          const Real min, const Real max) const;

    std::pair<std::vector<domain_id_type>, std::pair<domain_id_type, Real> >
    get_intrusive_domains(const vertex_id_type& vid,
                          const Real min, const Real max) const;

    //! calculate max circle size allowed by geometric restraints
    Real get_max_circle_size(const std::pair<Real3, face_id_type>& pos) const;
    //! calculate max cone size allowed by geometric restraints
    Real get_max_cone_size(const vertex_id_type& vid) const;

//     domain_type create_domain(
//             const ParticleID& pid, const Particle& p, const face_id_type fid);
//
//     Single create_single(
//             const ParticleID& pid, const Particle& p, const face_id_type fid);
//
//     circular_shell_type draw_circular_shell(
//             const ParticleID& pid, const Particle& p, const face_id_type fid,
//             const std::pair<vertex_id_type, Real>& closest_vertex,
//             const std::pair<domain_id_type, Real>& closest_domain);
//
//     circular_shell_type draw_conical_shell(
//             const ParticleID& pid, const Particle& p, const face_id_type fid,
//             const std::pair<vertex_id_type, Real>& closest_vertex,
//             const std::pair<domain_id_type, Real>& closest_domain);
//
//     std::pair<bool, domain_type>
//     form_pair_or_multi(const ParticleID& pid,
//                        const std::vector<domain_type>& bursted)
//     {
//         // not implemented yet.
//         return std::make_pair(false, domain_type());
//     }

  private:

    static const Real single_circular_shell_factor;

  private:

    // from SimulatorBase
    // boost::shared_ptr<model_type> model_;
    // boost::shared_ptr<world_type> world_;
    // Integer num_steps_;
    Real dt_;
    Real bd_dt_factor_;
    scheduler_type                     scheduler_;
    shell_id_generator_type            shell_id_gen;
    shell_container_type               shell_container_;
    mutable_shell_visitor_applier_type mut_sh_vis_applier;
//     std::vector<std::pair<ReactionRule, reaction_info_type> > last_reactions_;
};

template<typename T>
const Real SGFRDSimulator<T>::single_circular_shell_factor = 1.1;

template<typename T>
void SGFRDSimulator<T>::step()
{
    this->set_time(this->scheduler_.next_time());
    fire_event(this->scheduler_.pop());
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
    for(std::vector<std::pair<ParticleID, Particle> >::const_iterator
        iter = ps.begin(); iter != ps.end(); ++iter)
    {
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
                sim.escape(dom);
                return;
            }
            case Single::REACTION:
            {
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
    domain_firer firer(*this);
    boost::apply_visitor(firer, ev.second->domain());

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
        return;
    }

    void operator()(const Pair& dom)
    {
        pair_shell_burster burster(sim, dom, result);
        boost::apply_visitor(burster, sim.get_shell(dom.shell_id()));
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
    // when it called, there are no intruders in the domain region.

    single_shell_escapement(SGFRDSimulator& s, const Single& d)
        : sim(s), dom(d)
    {}

    void operator()(const circular_shell_type& sh)
    {
        if(sh.size() == dom.particle().radius()) return;

        Particle   p   = dom.particle();
        ParticleID pid = dom.particle_id();

        const Real r     = sh.size() - p.radius();
        const Real theta = sim.uniform_real() * 2 * M_PI;
        const face_id_type   fid  = sim.get_face_id(pid);
        const triangle_type& face = sim.polygon().triangle_at(fid);
        const Real3 direction = rotate(theta, face.normal(), face.represent());

        std::pair<std::pair<Real3, face_id_type>, Real3> state =
            std::make_pair(/*position = */std::make_pair(p.position(), fid),
                       /*displacement = */direction * r / length(direction));

        unsigned int continue_count = 2;
        while(continue_count > 0)
        {
            state = sim.polygon().move_next_face(state.first, state.second);
            const Real3& disp = state.second;
            if(disp[0] == 0. && disp[1] == 0. && disp[2] == 0.) break;
            --continue_count;
        }
        if(continue_count == 0)
            std::cerr << "[WARNING] moving on face: precision lost" << std::endl;

        p.position() = state.first.first;
        sim.update_particle(pid, p, state.first.second);
        sim.create_event(pid, p, state.first.second);
        return ;
    }

    void operator()(const conical_surface_shell_type& sh)
    {
        Particle           p   = dom.particle();
        const ParticleID   pid = dom.particle_id();
        const face_id_type fid = sim.get_face_id(pid);

        const Real r     = sh.size() - p.radius();
        greens_functions::GreensFunction2DRefWedgeAbs
            gf(p.D(), length(p.position() - sh.position()),
               r,     sh.shape().apex_angle());
        const Real theta = gf.drawTheta(sim.uniform_real(), r, dom.dt());

        const std::pair<Real3, face_id_type> state =
            sim.polygon().rotate_around_vertex(std::make_pair(p.position(), fid),
                                               sh.structure_id(), r, theta);
        p.position() = state.first;
        sim.update_particle(pid, p, state.second);
        sim.create_event(pid, p, state.second);
        return;
    }

  private:

    SGFRDSimulator<T>& sim;
    Single const& dom;
};

template<typename T>
void SGFRDSimulator<T>::escape(const Single& domain)
{
    single_shell_escapement escapement(*this, domain);
    boost::apply_visitor(escapement, get_shell(domain.shell_id()));
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
    const std::pair<Real3, face_id_type> pos = std::make_pair(p.position(), fid);
    const Real min_circle_size = p.radius() * single_circular_shell_factor;

    std::pair<std::vector<vertex_id_type>, std::pair<vertex_id_type, Real>
        > const intrusive_vertices(get_intrusive_vertices(pos, min_circle_size));

    if(intrusive_vertices.first.empty())
    {
        const Real max_circle_size = get_max_circle_size(pos);
        std::pair<std::vector<domain_id_type>, std::pair<domain_id_type, Real>
            > const intrusive_domains(get_intrusive_domains(
                        pos, min_circle_size, max_circle_size));

        if(intrusive_domains.first.empty())
        {
            return add_event(create_single(create_single_circular_shell(
                        pos, intrusive_domains.second.second), pid, p));
        }
        else
        {// TODO burst and form pair or multi if needed
            std::cerr << "[WARNING] intrusive domains exist." << std::endl;
            return add_event(create_single(create_single_circular_shell(
                        pos, intrusive_domains.second.second), pid, p));
        }
    }
    else // conical surface shell
    {
        const vertex_id_type& vid = intrusive_vertices.second.first;
        const Real max_cone_size = get_max_cone_size(vid);
        std::pair<std::vector<domain_id_type>, std::pair<domain_id_type, Real>
            > const intrusive_domains(get_intrusive_domains(
                        vid, min_circle_size, max_cone_size));

        if(intrusive_domains.first.empty())
        {
            return add_event(create_single(create_single_conical_surface_shell(
                std::make_pair(pos.first, vid), intrusive_domains.second.second),
                        pid, p));
        }
        else
        {// TODO burst and form pair or multi if needed
            std::cerr << "[WARNING] intrusive domains exist." << std::endl;
            return add_event(create_single(create_single_conical_surface_shell(
                std::make_pair(pos.first, vid), intrusive_domains.second.second),
                        pid, p));
        }
    }
}

template<typename T>
std::pair<ShellID, typename SGFRDSimulator<T>::circle_type>
SGFRDSimulator<T>::create_single_circular_shell(
            const std::pair<Real3, face_id_type>& pos, const Real size)
{
    const ShellID id(shell_id_gen());
    const circle_type shape(size, pos.first,
                            polygon().triangle_at(pos.second).normal());
    shell_container_.add_shell(id, circular_shell_type(shape, pos.second), pos.second);
    return std::make_pair(id, shape);
}

template<typename T>
std::pair<ShellID, typename SGFRDSimulator<T>::conical_surface_type>
SGFRDSimulator<T>::create_single_conical_surface_shell(
            const std::pair<Real3, vertex_id_type>& pos, const Real size)
{
    const ShellID id(shell_id_gen());
    const conical_surface_type shape(
            polygon().vertex_at(pos.second).position,
            polygon().apex_angle(pos.second), size);

    shell_container_.add_shell(id, conical_surface_shell_type(shape, pos.second),
                               pos.second);
    return std::make_pair(id, shape);
}


template<typename T>
Single SGFRDSimulator<T>::create_single(
        const std::pair<ShellID, circle_type>& sh,
        const ParticleID& pid, const Particle& p)
{
    throw ecell4::NotImplemented("get_intrusive_vertices");
}

template<typename T>
Single SGFRDSimulator<T>::create_single(
        const std::pair<ShellID, conical_surface_type>& sh,
        const ParticleID& pid, const Particle& p)
{
    throw ecell4::NotImplemented("get_intrusive_vertices");
}

template<typename T>
std::pair<std::vector<typename SGFRDSimulator<T>::vertex_id_type>,
          std::pair<typename SGFRDSimulator<T>::vertex_id_type, Real> >
SGFRDSimulator<T>::get_intrusive_vertices(
        const std::pair<Real3, face_id_type>& pos, const Real radius) const
{
    throw ecell4::NotImplemented("get_intrusive_vertices");
}

template<typename T>
std::pair<std::vector<typename SGFRDSimulator<T>::domain_id_type>,
          std::pair<typename SGFRDSimulator<T>::domain_id_type, Real> >
SGFRDSimulator<T>::get_intrusive_domains(
        const std::pair<Real3, face_id_type>& pos, const Real min, const Real max) const
{
    throw ecell4::NotImplemented("get_intrusive_vertices");
}

template<typename T>
std::pair<std::vector<typename SGFRDSimulator<T>::domain_id_type>,
          std::pair<typename SGFRDSimulator<T>::domain_id_type, Real> >
SGFRDSimulator<T>::get_intrusive_domains(const vertex_id_type& vid,
                      const Real min, const Real max) const
{
    throw ecell4::NotImplemented("get_intrusive_vertices");
}

template<typename T>
Real SGFRDSimulator<T>::get_max_circle_size(
        const std::pair<Real3, face_id_type>& pos) const
{
    throw ecell4::NotImplemented("get_max_circle_size");
}

template<typename T>
Real SGFRDSimulator<T>::get_max_cone_size(const vertex_id_type& vid) const
{
    throw ecell4::NotImplemented("get_max_cone_size");
}






} // sgfrd
} // ecell4
#endif // ECELL4_SGFRD_SIMULATOR