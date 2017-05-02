#ifndef ECELL4_SGFRD_SIMULATOR
#define ECELL4_SGFRD_SIMULATOR

#include <ecell4/core/SimulatorBase.hpp>
#include "ShellContainer.hpp"
#include "ShellVisitorApplier.hpp"
#include "ShellVisitors.hpp"
#include "SGFRDEvent.hpp"
#include "SGFRDWorld.hpp"

namespace ecell4
{
namespace sgfrd
{

template<typename T_polygon_traits>
class SGFRDSimulator :
    public ecell4::SimulatorBase<ecell4::Model, SGFRDWorld<T_polygon_traits> >
{
  public:
    typedef T_polygon_traits    polygon_traits_type;
    typedef SGFRDEvent          event_type;
    typedef SGFRDEventScheduler scheduler_type;
    typedef DomainID            domain_id_type;
    typedef typename SGFRDEvent::domain_type domain_type; //actually a variant
    typedef ecell4::SimulatorBase<ecell4::Model, SGFRDWorld<polygon_traits_type>
            > base_type;
    typedef typename base_type::world_type world_type;
    typedef typename base_type::model_type model_type;
    typedef ShellContainer<polygon_traits_type> shell_container_type;
    typedef typename shell_container_type::shell_type shell_type; // variant
    typedef typename shell_container_type::shell_id_pair_type shell_id_pair_type;
    typedef typename shell_container_type::circular_shell_type circular_shell_type;
    typedef typename shell_container_type::conical_surface_shell_type
        conical_surface_shell_type;

    typedef mutable_shell_visitor_applier<polygon_traits_type>
            mutable_shell_visitor_applier_type;

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

  public:

    SGFRDSimulator(const boost::shared_ptr<world_type>& world,
                   const boost::shared_ptr<model_type>& model,
                   Real bd_dt_factor = 1e-5)
        : base_type(model, world), dt_(0), bd_dt_factor_(bd_dt_factor),
          shell_container_(world->polygon())
    {
        initialize();
    }

    SGFRDSimulator(boost::shared_ptr<world_type> world, Real bd_dt_factor = 1e-5)
        : base_type(world), dt_(0), bd_dt_factor_(bd_dt_factor),
          shell_container_(world->polygon())
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
//
//     boost::shared_ptr<RandomNumberGenerator> rng() {return this->world_->rng();}

  private:

    // make event from domain and push it into scheduler
    template<typename domainT>
    void add_event(const domainT& dom);

    // pop event from scheduler then fire the event
    void fire_event(const event_type& event);
    void fire_domain(const Single&);
    void fire_domain(const Pair&);
    void fire_domain(const Multi&);

    // single event
    void escape(const Single& domain);
    void reaction(const Single& domain);

    // propagate particle
    void propagate(const circular_shell_type& sh,
                   const Real r, const Real theta);
    void propagate(const conical_surface_shell_type& sh,
                   const Real r, const Real theta);

    // burst domain(related to an event)
    void burst_event(const event_type& event,
                     std::vector<domain_type>& result);
    void burst_non_multis(std::vector<domain_type>& domains);

    domain_type burst_domain(const Single& dom);
    domain_type burst_domain(const Pair&   dom);
    domain_type burst_domain(const Multi&  dom);

    std::pair<std::vector<vertex_id_type>, std::pair<vertex_id_type, Real> >
    get_intrusive_vertices(const std::pair<Real3, face_id_type>& pos,
                           const Real radius) const;

    std::pair<std::vector<domain_id_type>, std::pair<domain_id_type, Real> >
    get_intrusive_domains(const Real3& pos, const Real radius) const;

    Single create_closely_fitted_single_circular_domain(
            const ParticleID& pid, const Particle& p, const face_id_type fid);

    domain_type create_domain(
            const ParticleID& pid, const Particle& p, const face_id_type fid);

    Single create_single(
            const ParticleID& pid, const Particle& p, const face_id_type fid);

    circular_shell_type draw_circular_shell(
            const ParticleID& pid, const Particle& p, const face_id_type fid,
            const std::pair<vertex_id_type, Real>& closest_vertex,
            const std::pair<domain_id_type, Real>& closest_domain);

    circular_shell_type draw_conical_shell(
            const ParticleID& pid, const Particle& p, const face_id_type fid,
            const std::pair<vertex_id_type, Real>& closest_vertex,
            const std::pair<domain_id_type, Real>& closest_domain);

    std::pair<bool, domain_type>
    form_pair_or_multi(const ParticleID& pid,
                       const std::vector<domain_type>& bursted)
    {
        // not implemented yet.
        return std::make_pair(false, domain_type());
    }

  private:

    // consider scheduler.time
    Real time() const {return this->world_->t();}
    void set_time(const Real t) {return this->world_->set_t(t);}

  private:

    // from SimulatorBase
    // boost::shared_ptr<model_type> model_;
    // boost::shared_ptr<world_type> world_;
    // Integer num_steps_;
    Real dt_;
    Real bd_dt_factor_;
    scheduler_type       scheduler_;
    shell_container_type shell_container_;
//     std::vector<std::pair<ReactionRule, reaction_info_type> > last_reactions_;
};

template<typename T>
void SGFRDSimulator<T>::step()
{
    this->set_time(this->scheduler_.next_time());
    fire_event(*(this->scheduler_.pop().second));
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
    //XXX copy of all the particles occur
    std::vector<std::pair<ParticleID, Particle> > ps =
        this->world_->list_particles();
    for(std::vector<std::pair<ParticleID, Particle> >::const_iterator
        iter = ps.begin(); iter != ps.end(); ++iter)
    {
        this->add_event(this->create_closely_fitted_single_circular_domain(
                    iter->first, iter->second,
                    this->world_->get_faceID(iter->first)));
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
    const DomainID did =
        scheduler_.add(boost::make_shared<event_type>(this->time(), dom));
    mutable_shell_visitor_applier_type applier(this->shell_container_);
    applier(domain_id_setter(did), dom);
    return ;
}





} // sgfrd
} // ecell4
#endif // ECELL4_SGFRD_SIMULATOR
