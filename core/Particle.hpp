#ifndef __ECELL4_PARTICLE_HPP
#define __ECELL4_PARTICLE_HPP

#include <map>

#include "types.hpp"
#include "Position3.hpp"
#include "Species.hpp"
#include "Identifier.hpp"
#include "config.h"


namespace ecell4
{

class Particle
{
public:
    typedef Position3 position_type;
    typedef Real length_type;
    typedef Real D_type;
    typedef Species species_type;

public:
    Particle()
    {
        ;
    }

    Particle(
        const Species& sp, const Position3& pos, const Real& radius,
        const Real& D)
        : species_(sp), position_(pos), radius_(radius), D_(D)
    {
        ;
    }

    Position3& position()
    {
        return position_;
    }

    const Position3& position() const
    {
        return position_;
    }

    Real& radius()
    {
        return radius_;
    }

    const Real& radius() const
    {
        return radius_;
    }

    Real& D()
    {
        return D_;
    }

    const Real& D() const
    {
        return D_;
    }

    Species& species()
    {
        return species_;
    }

    /*
    Species::serial_type &sid()
    {
        return species_.serial();
    }*/

    const Species::serial_type &sid() const
    {
        return species_.serial();
    }

    const Species& species() const
    {
        return species_;
    }

private:

    Species species_;
    Position3 position_;
    Real radius_, D_;
};

template<typename Tstrm_, typename Ttraits_>
inline std::basic_ostream<Tstrm_, Ttraits_>& operator<<(std::basic_ostream<Tstrm_, Ttraits_>& strm, const Particle& p)
{
    strm << "Particle(" << "{ " << p.position() << ", " << p.radius() << "}, " << ", D=" << p.D() << ", " << p.sid() << ")";
    return strm;
}

} // ecell4


#if defined(HAVE_TR1_FUNCTIONAL)
namespace std { namespace tr1 {
#elif defined(HAVE_STD_HASH)
namespace std {
#elif defined(HAVE_BOOST_FUNCTIONAL_HASH_HPP)
namespace boost {
#endif

template<>
struct hash<ecell4::Particle>
{
    typedef ecell4::Particle argument_type;

    std::size_t operator()(argument_type const& val)
    {
        return hash<typename argument_type::position_type>()(val.position()) ^
            hash<argument_type::length_type>()(val.radius()) ^
            hash<argument_type::D_type>()(val.D()) ^
            hash<argument_type::species_type>()(val.species());
    }
};

#if defined(HAVE_TR1_FUNCTIONAL)
} } // namespace std::tr1
#elif defined(HAVE_STD_HASH)
} // namespace std
#elif defined(HAVE_BOOST_FUNCTIONAL_HASH_HPP)
} // namespace boost
#endif

#endif /* __ECELL4_PARTICLE_HPP */
