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
template<typename T_polygon_traits>
class Single
{
  public:

    typedef T_polygon_traits polygon_traits;
    typedef typename polygon_traits::face_id_type face_id_type;
    typedef typename polygon_traits::vertex_id_type vertex_id_type;

    typedef DomainID   identifier_type;
    typedef Particle   particle_type;
    typedef ParticleID particle_id_type;
    typedef std::pair<ParticleID, Particle> particle_id_pair;
    typedef Shell<ecell4::Circle, face_id_type> circular_shell;
    typedef Shell<ecell4::ConicalSurface, vertex_id_type> conical_surface_shell;
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
    identifier_type    id_;
    particle_id_pair   particle_;
    storage_type       shell_;
};


} // sgfrd
} // ecell4
#endif /* ECELL4_SGFRD_SINGLE_DOMAIN */
