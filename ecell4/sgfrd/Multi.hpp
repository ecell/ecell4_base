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

    Multi(simulator_type& sim, world_type& world,
          Real dt_factor = 0.01, Real rl_factor = 0.1 /* [0.05 ~ 0.1] */)
        : kind_(NONE), dt_(-1.0), begin_time_(0.), reaction_length_(-1.0),
          dt_factor_(dt_factor), reaction_length_factor_(rl_factor),
          simulator_(sim), world_(world), model_(*world.lock_model()),
          container_(world)
    {}

    Multi(simulator_type& sim, world_type& world, Real begin_t,
          Real dt_factor = 0.01, Real rl_factor = 0.1 /* [0.05 ~ 0.1] */)
        : kind_(NONE), dt_(-1.0), begin_time_(begin_t), reaction_length_(-1.0),
          dt_factor_(dt_factor), reaction_length_factor_(rl_factor),
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
        assert(this->dt_              > 0.0);
        assert(this->reaction_length_ > 0.0);

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

    void determine_parameters()
    {
        this->determine_reaction_length();
        this->determine_delta_t();
        return;
    }

    void determine_reaction_length()
    {
        Real r_min =  std::numeric_limits<Real>::max();
        for(const auto& p : this->particles())
        {
            r_min = std::min(p.second.radius(), r_min);
        }
        reaction_length_ = this->reaction_length_factor_ * r_min;
        return;
    }
    void determine_delta_t()
    {
        // it assumes this->determine_reaction_length() has already been called
        assert(this->reaction_length_ > 0.0);

        // collect possible reactions and find maximum D
        std::vector<Species> sps;
        sps.reserve(this->particles().size());
        Real D_max = -std::numeric_limits<Real>::max();
        for(const auto& p : this->particles())
        {
            D_max = std::max(p.second.D(), D_max);
            if(std::find(sps.begin(), sps.end(), p.second.species()) == sps.end())
            {
                sps.push_back(p.second.species());
            }
        }

        Real k_max = 0.0;
        if(sps.size() >= 2)
        {
            // TODO: make it more efficient ... ?
            for(std::size_t i=0; i<sps.size(); ++i)
            {
                const auto& sp1 = sps.at(i);
                // start from j = i. there can be a reaction like "A+A -> B".
                for(std::size_t j=i; j<sps.size(); ++j)
                {
                    const auto& sp2  = sps.at(j);
                    for(auto&& rule : model_.query_reaction_rules(sp1, sp2))
                    {
                        k_max = std::max(k_max, rule.k());
                    }
                }
            }
        }

        const Real P_max = 0.01; // upper limit for reaction probability
        const Real delta = this->reaction_length_;

        const Real upper_limit_D = delta * delta / D_max;
        const Real upper_limit_k = (k_max > 0.0) ? (P_max * delta / k_max) :
                                   std::numeric_limits<Real>::infinity();
        this->dt_ = this->dt_factor_ * std::min(upper_limit_D, upper_limit_k);
        return;
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
    Real dt_factor_, reaction_length_factor_;
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
