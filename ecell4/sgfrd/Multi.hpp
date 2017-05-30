#ifndef ECELL4_SGFRD_MULTI_DOMAIN
#define ECELL4_SGFRD_MULTI_DOMAIN
#include "SGFRDWorld.hpp"
#include <ecell4/sgfrd/ShellID.hpp>
#include <ecell4/sgfrd/BDPropagator.hpp>
#include <ecell4/core/Particle.hpp>

namespace ecell4
{
namespace sgfrd
{

class SGFRDSimulator;
class SGFRDWorld;

class Multi
{
  public:

    enum EventKind
    {
        NONE,
        ESCAPE,
        REACTION,
    };

    typedef polygon_traits polygon_traits_type;
    typedef ecell4::Polygon<polygon_traits_type>     polygon_type;

    typedef Particle   particle_type;
    typedef ParticleID particle_id_type;
    typedef ShellID    shell_id_type;
    typedef std::vector<particle_id_type> particle_ids_type;
    typedef std::vector<shell_id_type>    shell_ids_type;
    // TODO implement multi-particle-container
    typedef BDPropagator<SGFRDWorld>      propagator_type;
    typedef SGFRDWorld                    world_type;
    typedef SGFRDSimulator                simulator_type;

  public:
    Multi(simulator_type& sim, world_type& world)
        : dt_(0.), begin_time_(0.), simulator_(sim), world_(world) {}
    ~Multi(){}

    void step(){/*TODO*/ return;}

    EventKind& eventkind()       {return kind_;}
    EventKind  eventkind() const {return kind_;}

    Real& dt()       {return dt_;}
    Real  dt() const {return dt_;}
    Real& begin_time()       {return begin_time_;}
    Real  begin_time() const {return begin_time_;}

    void add_particle(particle_id_type const& pid){particles_.push_back(pid);}
    void add_shell   (shell_id_type    const& sid){shells_.push_back(sid);}

    shell_ids_type&          shell_ids()       {return shells_;}
    shell_ids_type const&    shell_ids() const {return shells_;}
    particle_ids_type&       particle_ids()       {return particles_;}
    particle_ids_type const& particle_ids() const {return particles_;}

    std::size_t num_shells()   const {return shells_.size();}
    std::size_t multiplicity() const {return particles_.size();}

  private:

    EventKind kind_;
    Real dt_, begin_time_;
    simulator_type&   simulator_;
    world_type&       world_;
    particle_ids_type particles_;
    shell_ids_type    shells_;
};

} // sgfrd
} // ecell4
#endif /* ECELL4_SGFRD_MULTI_DOMAIN */
