#ifndef ECELL4_SGFRD_SINGLE_DOMAIN
#define ECELL4_SGFRD_SINGLE_DOMAIN
#include <ecell4/sgfrd/ShellID.hpp>
#include <ecell4/core/Particle.hpp>

namespace ecell4
{
namespace sgfrd
{

class Single
{
  public:

    enum EventKind
    {
        ESCAPE,
        REACTION,
        UNKNOWN,
    };

    typedef ShellID    shell_id_type;
    typedef Particle   particle_type;
    typedef ParticleID particle_id_type;
    typedef std::pair<ParticleID, Particle> particle_id_pair_type;

  public:
    Single(): kind_(UNKNOWN), dt_(0.), begin_time_(0.){}
    Single(const EventKind kind, const Real dt, const Real begin_time,
           const shell_id_type shid, const particle_id_pair_type& pidp)
        : kind_(kind), dt_(dt), begin_time_(begin_time),
          shell_id_(shid), particle_(pidp)
    {}
    ~Single(){}

    shell_id_type&       shell_id()       {return shell_id_;}
    shell_id_type const& shell_id() const {return shell_id_;}

    Real& dt()       {return dt_;}
    Real  dt() const {return dt_;}
    Real& begin_time()       {return begin_time_;}
    Real  begin_time() const {return begin_time_;}

    EventKind  eventkind() const {return kind_;}
    EventKind& eventkind()       {return kind_;}

    particle_id_pair_type&       particle_id_pair()       {return particle_;}
    particle_id_pair_type const& particle_id_pair() const {return particle_;}

    Particle&       particle()       {return particle_.second;}
    Particle const& particle() const {return particle_.second;}
    ParticleID&       particle_id()       {return particle_.first;}
    ParticleID const& particle_id() const {return particle_.first;}

    std::size_t num_shells() const {return 1;}
    std::size_t multiplicity() const {return 1;}

  private:

    EventKind kind_;
    Real dt_;
    Real begin_time_;
    shell_id_type    shell_id_;
    particle_id_pair_type particle_;
};


} // sgfrd
} // ecell4
#endif /* ECELL4_SGFRD_SINGLE_DOMAIN */
