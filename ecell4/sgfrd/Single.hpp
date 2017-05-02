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

    typedef ShellID    shell_id_type;
    typedef Particle   particle_type;
    typedef ParticleID particle_id_type;
    typedef std::pair<ParticleID, Particle> particle_id_pair;

  public:
    Single(): dt_(0.), last_time_(0.){}
    Single(const Real dt, const Real last_time, const shell_id_type& sh)
        : dt_(dt), last_time_(last_time), shell_id_(sh)
    {}
    ~Single(){}

    shell_id_type&       shell_id()       {return shell_id_;}
    shell_id_type const& shell_id() const {return shell_id_;}

    Real& dt()       {return dt_;}
    Real  dt() const {return dt_;}
    Real& last_time()       {return last_time_;}
    Real  last_time() const {return last_time_;}

    particle_id_pair&       particle()       {return particle_;}
    particle_id_pair const& particle() const {return particle_;}

    std::size_t num_shells() const {return 1;}
    std::size_t multiplicity() const {return 1;}

  private:

    Real dt_;
    Real last_time_;
    shell_id_type    shell_id_;
    particle_id_pair particle_;
};


} // sgfrd
} // ecell4
#endif /* ECELL4_SGFRD_SINGLE_DOMAIN */
