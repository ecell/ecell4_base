#ifndef __ECELL4_PARTICLE_HPP
#define __ECELL4_PARTICLE_HPP

#include "config.h"

#include <ostream>
#if defined(HAVE_TR1_FUNCTIONAL)
#include <tr1/functional>
#elif defined(HAVE_STD_HASH)
#include <functional>
#elif defined(HAVE_BOOST_FUNCTIONAL_HASH_HPP)
#include <boost/functional/hash.hpp>
#endif

#include <map>

#include "types.hpp"
#include "Position3.hpp"
#include "Species.hpp"
#include "Identifier.hpp"


namespace ecell4
{

class Particle
{
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

    const Species& species() const
    {
        return species_;
    }

private:

    Species species_;
    Position3 position_;
    Real radius_, D_;
};

struct ParticleID:
        public Identifier<ParticleID, unsigned long long, int>
{
    typedef Identifier<ParticleID, unsigned long long, int> base_type;

    ParticleID(const value_type& value = value_type(0, 0))
        : base_type(value)
    {
        ;
    }
};

} // ecell4

#if defined(HAVE_TR1_FUNCTIONAL)
namespace std
{

namespace tr1
{
#elif defined(HAVE_STD_HASH)
namespace std
{
#elif defined(HAVE_BOOST_FUNCTIONAL_HASH_HPP)
namespace boost
{
#endif

template<>
struct hash<ecell4::ParticleID>
{
    std::size_t operator()(const ecell4::ParticleID& val) const
    {
        return static_cast<std::size_t>(val().first ^ val().second);
    }
};

#if defined(HAVE_TR1_FUNCTIONAL)
} // tr1

} // std
#elif defined(HAVE_STD_HASH)
} // std
#elif defined(HAVE_BOOST_FUNCTIONAL_HASH_HPP)
} // boost
#endif

#endif /* __ECELL4_PARTICLE_HPP */
