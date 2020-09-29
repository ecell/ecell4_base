#ifndef ECELL4_NGFRD_SINGLE_DOMAIN
#define ECELL4_NGFRD_SINGLE_DOMAIN
#include <ecell4/ngfrd/ShellID.hpp>
#include <ecell4/core/Identifier.hpp>
#include <cstdint>

namespace ecell4
{
namespace ngfrd
{

class SingleDomain
{
  public:

    enum class EventKind : std::uint8_t
    {
        Escape,
        Reaction,
        Unknown,
    };
    using shell_id_type    = ShellID;
    using particle_id_type = ParticleID;

  public:

    SingleDomain(): kind_(EventKind::Unknown), dt_(0.), begin_time_(0.){}
    SingleDomain(const EventKind kind, const Real dt, const Real begin_time,
           const shell_id_type shid, const particle_id_type& pid)
        : kind_(kind), dt_(dt), begin_time_(begin_time),
          shell_id_(shid), particle_id_(pid)
    {}
    ~SingleDomain() = default;

    shell_id_type&       shell_id()       noexcept {return shell_id_;}
    shell_id_type const& shell_id() const noexcept {return shell_id_;}
    ParticleID&       particle_id()       noexcept {return particle_id_;}
    ParticleID const& particle_id() const noexcept {return particle_id_;}

    Real& dt()       noexcept {return dt_;}
    Real  dt() const noexcept {return dt_;}
    Real& begin_time()       noexcept {return begin_time_;}
    Real  begin_time() const noexcept {return begin_time_;}

    EventKind  eventkind() const {return kind_;}
    EventKind& eventkind()       {return kind_;}

    constexpr std::size_t num_shells()   const noexcept {return 1;}
    constexpr std::size_t multiplicity() const noexcept {return 1;}

  private:

    EventKind kind_;
    Real dt_;
    Real begin_time_;
    shell_id_type    shell_id_;
    particle_id_type particle_id_;
};

} // sgfrd
} // ecell4
#endif /* ECELL4_SGFRD_SINGLE_DOMAIN */
