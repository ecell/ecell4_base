#ifndef ECELL4_SGFRD_SINGLE_DOMAIN
#define ECELL4_SGFRD_SINGLE_DOMAIN
#include <ecell4/sgfrd/ShellID.hpp>
#include <ecell4/sgfrd/DomainID.hpp>
#include <ecell4/core/Particle.hpp>

namespace ecell4
{
namespace sgfrd
{

class Single
{
  public:

    typedef T_polygon_traits polygon_traits;
    typedef typename polygon_traits::face_id_type face_id_type;
    typedef typename polygon_traits::vertex_id_type vertex_id_type;

    typedef ShellID    shell_id_type;
    typedef DomainID   identifier_type;
    typedef Particle   particle_type;
    typedef ParticleID particle_id_type;
    typedef std::pair<ParticleID, Particle> particle_id_pair;

  public:
    Single(): dt_(0.), last_time_(0.){}
    Single(const identifier_type& id, const shell_id_type& sh)
        : id_(id), shell_id_(sh)
    {}
    ~Single(){}

    identifier_type&       id()       {return id_;}
    identifier_type const& id() const {return id_;}
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
    identifier_type  id_;
    shell_id_type    shell_id_;
    particle_id_pair particle_;
};


} // sgfrd
} // ecell4
#endif /* ECELL4_SGFRD_SINGLE_DOMAIN */
