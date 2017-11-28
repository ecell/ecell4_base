#ifndef ECELL4_SGFRD_PAIR_DOMAIN
#define ECELL4_SGFRD_PAIR_DOMAIN
#include <ecell4/sgfrd/ShellID.hpp>
#include <ecell4/core/Particle.hpp>

namespace ecell4
{
namespace sgfrd
{

class Pair
{
  public:

    enum EventKind
    {
        SINGLE_REACTION_1,
        SINGLE_REACTION_2,
        COM_ESCAPE,
        IV_ESCAPE,
        IV_REACTION,
        IV_UNDETERMINED,
        UNDETERMINED
    };

    typedef ShellID    shell_id_type;
    typedef Particle   particle_type;
    typedef ParticleID particle_id_type;
    typedef std::pair<ParticleID, Particle>   particle_id_pair;
    typedef boost::array<particle_id_pair, 2> particle_array_type;

    static Real calc_D_ipv(const Real D1, const Real D2) throw()
    {
        return D1 + D2;
    }

    static Real calc_D_com(const Real D1, const Real D2) throw()
    {
        return D1 * D2 / (D1 + D2);
    }

    // TODO: calculate more efficient R_*
    static Real calc_R_ipv(const Real r_shell,
                           const Particle& p1, const Particle& p2) throw()
    {
        return r_shell / 2;
    }
    static Real calc_R_com(const Real r_shell,
                           const Particle& p1, const Particle& p2) throw()
    {
        return r_shell / 2;
    }

  public:

    Pair(): kind_(UNDETERMINED), dt_(0.), begin_time_(0.){}
    ~Pair(){}

    Pair(const EventKind kind, const Real dt, const Real begin_time,
         const shell_id_type& sh, const Real shell_rad,
         const particle_id_pair& p0, const particle_id_pair& p1)
        : kind_(kind), dt_(dt), begin_time_(begin_time), shell_id_(sh)
    {
        particles_[0] = p0;
        particles_[1] = p1;
        this->r_ipv_ = Pair::calc_R_ipv(shell_rad, p0.second, p1.second);
        this->r_com_ = Pair::calc_R_com(shell_rad, p0.second, p1.second);
    }

    shell_id_type&       shell_id()       throw() {return shell_id_;}
    shell_id_type const& shell_id() const throw() {return shell_id_;}

    Real& dt()               throw() {return dt_;}
    Real  dt()         const throw() {return dt_;}
    Real& begin_time()       throw() {return begin_time_;}
    Real  begin_time() const throw() {return begin_time_;}

    Real R_ipv() const throw() {return this->r_ipv_;}
    Real R_com() const throw() {return this->r_com_;}

    particle_id_pair&       operator[](std::size_t i)       throw()
    {return particles_[i];}
    particle_id_pair const& operator[](std::size_t i) const throw()
    {return particles_[i];}

    std::size_t num_shells()   const throw() {return 1;}
    std::size_t multiplicity() const throw() {return 2;}

  private:

    EventKind kind_;
    Real      dt_;
    Real      begin_time_;
    Real      r_ipv_;
    Real      r_com_;
    shell_id_type       shell_id_;
    particle_array_type particles_;
};

} // sgfrd
} // ecell4
#endif /* ECELL4_SGFRD_PAIR_DOMAIN */
