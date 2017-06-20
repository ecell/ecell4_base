#ifndef ECELL4_SGFRD_MULTI_DOMAIN
#define ECELL4_SGFRD_MULTI_DOMAIN
#include <ecell4/sgfrd/ShellID.hpp>
#include <ecell4/sgfrd/BDPropagator.hpp>
#include <ecell4/sgfrd/Informations.hpp>
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

    typedef Particle   particle_type;
    typedef ParticleID particle_id_type;
    typedef ShellID    shell_id_type;
    typedef std::pair<ParticleID, Particle>    particle_id_pair_type;
    typedef std::vector<particle_id_pair_type> particles_type;
    typedef std::vector<shell_id_type>         shell_ids_type;

    typedef SGFRDWorld               world_type;
    typedef world_type::polygon_type polygon_type;
    typedef world_type::model_type   model_type;
    typedef SGFRDSimulator           simulator_type;
    typedef MultiContainer           container_type;

    typedef ecell4::ReactionRule reaction_rule_type;
    typedef ecell4::sgfrd::MoleculeInfo molecule_info_type;
    typedef ecell4::sgfrd::ReactionInfo reaction_info_type;
    typedef std::pair<reaction_rule_type, reaction_info_type> reaction_log_type;
    typedef std::vector<reaction_log_type>                reaction_archive_type;

  public:

    Multi(simulator_type& sim, world_type& world)
        : dt_(1e-5), begin_time_(0.), reaction_length_(1e-3),
          simulator_(sim), world_(world),
          container_(world), model_(*world.lock_model())
    {}
    ~Multi(){}

    template<typename vcT>
    void step(vcT vc)
    {
        this->last_reactions_.clear();
        kind_ = NONE;

        BDPropagator<container_type, vcT> propagator(model_, container_,
                *(world_.polygon()), *(world_.rng()), dt_, reaction_length_,
                last_reactions_, vc);

        while(propagator())
        {
            // if reaction occurs, return immediately
            if(!last_reactions_.empty())
            {
                kind_ = REACTION;
                break;
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

    // TODO
    void determine_reaction_length(){return;}
    void determine_delta_t(){return;}

    void add_particle(particle_id_type const& pid){container_.make_entry(pid);}
    void add_shell   (shell_id_type    const& sid){shells_.push_back(sid);}

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
