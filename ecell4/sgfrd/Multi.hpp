#ifndef ECELL4_SGFRD_DOMAIN
#define ECELL4_SGFRD_DOMAIN
#include <ecell4/sgfrd/ShellID.hpp>
#include <ecell4/core/Particle.hpp>

namespace ecell4
{
namespace sgfrd
{

class Multi
{
  public:

//     typedef BDPropagator simulator_type;
    typedef Particle   particle_type;
    typedef ParticleID particle_id_type;
    typedef ShellID    shell_id_type;
    typedef std::pair<ParticleID, Particle> particle_id_pair;
    typedef std::vector<particle_id_pair>   particle_container_type;
    typedef std::vector<shell_id_type>      shell_id_container_type;

  public:
    Multi(): dt_(0.), begin_time_(0.){}
    ~Multi(){}

//     Multi(simulator_type& sim) : simulator_(sim){}

    Real& dt()       {return dt_;}
    Real  dt() const {return dt_;}
    Real& begin_time()       {return begin_time_;}
    Real  begin_time() const {return begin_time_;}

    shell_id_container_type&          shells()       {return shells_;}
    shell_id_container_type const&    shells() const {return shells_;}
    particle_container_type&       particles()       {return particles_;}
    particle_container_type const& particles() const {return particles_;}

    void add_particle(particle_id_pair const& pp)
    {
        particles_.push_back(pp);
        return;
    }
    void add_shell(shell_id_type const& sh)
    {
        shells_.push_back(sh);
        return;
    }

    std::size_t num_shells()   const {return shells_.size();}
    std::size_t multiplicity() const {return particles_.size();}

  private:

    Real dt_factor_;
    Real dt_, begin_time_;
//     simulator_type&         simulator_;
    particle_container_type particles_;
    shell_id_container_type    shells_;
};

} // sgfrd
} // ecell4
#endif /* ECELL4_SGFRD_MULTI_DOMAIN */
