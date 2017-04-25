#ifndef ECELL4_SGFRD_DOMAIN
#define ECELL4_SGFRD_DOMAIN
#include <ecell4/sgfrd/Shell.hpp>
#include <ecell4/core/Particle.hpp>
#include <ecell4/core/Identifier.hpp>
#include <ecell4/core/Circle.hpp>
#include <ecell4/core/Cone.hpp>
#include <boost/variant.hpp>

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
    typedef std::pair<ParticleID, Particle> particle_id_pair;
    typedef std::vector<particle_id_pair>   particle_container_type;

    typedef ShellID shell_id_type;
    typedef Shell<ecell4::Circle>           circular_shell;
    typedef Shell<ecell4::ConicalSurface>   conical_surface_shell;
    typedef boost::variant<circular_shell, conical_surface_shell> storage_type;
    typedef std::map<shell_id_type, storage_type> shell_container_type;

  public:
    Multi(){}
    ~Multi(){}

//     Multi(identifier_type const& id, simulator_type& sim)
//         : id_(id), simulator_(sim)
//     {}

    container_type&            shells()          {return shells_;}
    container_type const&      shells()    const {return shells_;}
    particle_array_type&       particles()       {return particle_;}
    particle_array_type const& particles() const {return particle_;}

    void add_particle(particle_id_pair const& pp)
    {
        particles_.update_particle(pp);
        return;
    }
    void add_shell(circular_shell_id_pair const& sh)
    {
        typename shell_container_type::value_type newsh(
                sh.first, storage_type(sh.second));
        boost::get<circular_shell>(newsh.second).domain_id() = this->id_;
        shells_.insert(newsh);
        return;
    }
    void add_shell(conical_shell_id_pair const& sh)
    {
        typename shell_container_type::value_type newsh(
                sh.first, storage_type(sh.second));
        boost::get<conical_shell>(newsh.second).domain_id() = this->id_;
        shells_.insert(newsh);
        return;
    }

    std::size_t num_shells()   const {return shells_.size();}
    std::size_t multiplicity() const {return particles_.size();}

  private:

    Real dt_factor_;
//     simulator_type&         simulator_;
    event_id_pair_type      event_;
    particle_container_type particles_;
    shell_container_type    shells_;
};

} // sgfrd
} // ecell4
#endif /* ECELL4_SGFRD_MULTI_DOMAIN */
