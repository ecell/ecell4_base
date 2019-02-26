#ifndef ECELL4_SGFRD_MULTI_DOMAIN
#define ECELL4_SGFRD_MULTI_DOMAIN
#include <ecell4/sgfrd/ShellID.hpp>
#include <ecell4/sgfrd/BDPropagator.hpp>
#include <ecell4/sgfrd/ReactionInfo.hpp>
#include <ecell4/sgfrd/MultiContainer.hpp>
#include <ecell4/core/Particle.hpp>

namespace ecell4
{
namespace sgfrd
{

class SGFRDSimulator;

class Multi
{
  public:

    static Real bd_factor() noexcept {return 0.05;}

    enum EventKind
    {
        NONE,
        ESCAPE,
        REACTION,
    };

    typedef std::pair<ParticleID, Particle>    particle_id_pair_type;
    typedef std::vector<particle_id_pair_type> particles_type;
    typedef std::vector<ShellID>               shell_ids_type;

    typedef SGFRDWorld               world_type;
    typedef world_type::polygon_type polygon_type;
    typedef world_type::model_type   model_type;
    typedef SGFRDSimulator           simulator_type;
    typedef MultiContainer           container_type;

    typedef ecell4::ReactionRule reaction_rule_type;
    typedef ecell4::sgfrd::ReactionInfo reaction_info_type;
    typedef std::pair<reaction_rule_type, reaction_info_type> reaction_log_type;
    typedef std::vector<reaction_log_type>                reaction_archive_type;

  public:

    Multi(simulator_type& sim, world_type& world, Real dt, Real rl)
        : kind_(NONE), dt_(dt), begin_time_(0.), reaction_length_(rl),
          simulator_(sim), world_(world), model_(*world.lock_model()),
          container_(world)
    {}

    Multi(simulator_type& sim, world_type& world, Real begin_t, Real dt, Real rl)
        : kind_(NONE), dt_(dt), begin_time_(begin_t), reaction_length_(rl),
          simulator_(sim), world_(world), model_(*world.lock_model()),
          container_(world)
    {}
    ~Multi(){}

    template<typename vcT>
    void step(vcT vc)
    {
        step(vc, this->dt_);
        return ;
    }

    template<typename vcT>
    void step(vcT vc, const Real dt)
    {
        this->last_reactions_.clear();
        kind_ = NONE;

        BDPropagator<container_type, vcT> propagator(model_, container_,
                *(world_.polygon()), *(world_.rng()), dt, reaction_length_,
                last_reactions_, vc);

        while(propagator())
        {
            // if reaction occurs, return immediately
            // XXX is it okay?
            if(!last_reactions_.empty())
            {
                kind_ = REACTION;
                break;
            }
            if(propagator.vc().escaped())
            {
                kind_ = ESCAPE;
            }
        }
    }

    EventKind& eventkind()       {return kind_;}
    EventKind  eventkind() const {return kind_;}

    Real& dt()       {return dt_;}
    Real  dt() const {return dt_;}
    Real& begin_time()       {return begin_time_;}
    Real  begin_time() const {return begin_time_;}
    Real& reaction_length()       {return reaction_length_;}
    Real  reaction_length() const {return reaction_length_;}

    void determine_reaction_length()
    {
//        Real r_min =  std::numeric_limits<Real>::max();
//        for(const auto& p : this->particles())
//        {
//            r_min = std::min(p.second.radius(), r_min);
//        }
//        reaction_length_ = bd_factor() * r_min;
//        return;
    }
    void determine_delta_t()
    {
//        Real r_min =  std::numeric_limits<Real>::max();
//        Real D_max = -std::numeric_limits<Real>::max();
//        for(const auto& p : this->particles())
//        {
//            r_min = std::min(p.second.radius(), r_min);
//            D_max = std::max(p.second.D(), D_max);
//        }
//        dt_ = (bd_factor() * r_min) * (bd_factor() * r_min) / D_max;
//        return;
    }

    bool add_particle(ParticleID const& pid)
    {
        return container_.make_entry(pid);
    }
    bool add_shell(ShellID const& sid)
    {
        if(std::find(shells_.begin(), shells_.end(), sid) != shells_.end())
        {
            return false;
        }
        shells_.push_back(sid);
        return true;
    }

    shell_ids_type&       shell_ids()       {return shells_;}
    shell_ids_type const& shell_ids() const {return shells_;}
    particles_type&       particles()       {return container_.list_particles();}
    particles_type const& particles() const {return container_.list_particles();}

    std::size_t num_shells()   const {return shells_.size();}
    std::size_t multiplicity() const {return container_.num_particles();}

    reaction_archive_type const& last_reactions() const {return last_reactions_;}

  private:

    EventKind kind_;
    Real dt_, begin_time_, reaction_length_;
    simulator_type&       simulator_;
    world_type&           world_;
    model_type&           model_;
    container_type        container_;
    shell_ids_type        shells_;
    reaction_archive_type last_reactions_;
};

} // sgfrd
} // ecell4
#endif /* ECELL4_SGFRD_MULTI_DOMAIN */
