#ifndef ECELL4_SGFRD_SINGLE_DOMAIN
#define ECELL4_SGFRD_SINGLE_DOMAIN
#include <ecell4/sgfrd/Shell.hpp>
#include <ecell4/sgfrd/DomainID.hpp>
#include <ecell4/core/Particle.hpp>
#include <ecell4/core/Identifier.hpp>
#include <ecell4/core/Circle.hpp>
#include <ecell4/core/Cone.hpp>
#include <boost/variant.hpp>

namespace ecell4
{
namespace sgfrd
{

/*! @brief single domain type. */
class Single
{
  public:

    typedef DomainID   identifier_type;
    typedef Particle   particle_type;
    typedef ParticleID particle_id_type;
    typedef std::pair<ParticleID, Particle> particle_id_pair;
    typedef Shell<ecell4::Circle>         circular_shell;
    typedef Shell<ecell4::ConicalSurface> conical_surface_shell;
    typedef boost::variant<circular_shell, conical_surface_shell> storage_type;

  public:
    Single(): dt_(0.), last_time_(0.){}
    ~Single(){}

    explicit Single(identifier_type const& id): id_(id){}

    Single(identifier_type const& id, circular_shell const& sh)
        : id_(id), shell_(sh)
    {}
    Single(identifier_type const& id, conical_surface_shell const& sh)
        : id_(id), shell_(sh)
    {}

    identifier_type&       id()       {return id_;}
    identifier_type const& id() const {return id_;}
    storage_type&       shell()       {return shell_;}
    storage_type const& shell() const {return shell_;}
    time_type& dt()       {return dt_;}
    time_type  dt() const {return dt_;}
    time_type& last_time()       {return last_time_;}
    time_type  last_time() const {return last_time_;}

    particle_id_pair&       particle()       {return particle_;}
    particle_id_pair const& particle() const {return particle_;}

    std::size_t num_shells() const {return 1;}
    std::size_t multiplicity() const {return 1;}

  private:

    Real dt_;
    Real last_time_;
    identifier_type    id_;
    event_id_pair_type event_;
    particle_id_pair   particle_;
    storage_type       shell_;
};


} // sgfrd
} // ecell4
#endif /* ECELL4_SGFRD_SINGLE_DOMAIN */
