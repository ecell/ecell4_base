#ifndef __PARTICLE_HPP
#define __PARTICLE_HPP

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
        Species const& sp, Position3 const& pos, Real const& radius,
        Real const& D)
        : species_(sp), position_(pos), radius_(radius), D_(D)
    {
        ;
    }

    Position3& position()
    {
        return position_;
    }

    Position3 const& position() const
    {
        return position_;
    }

    Real& radius()
    {
        return radius_;
    }

    Real const& radius() const
    {
        return radius_;
    }

    Real& D()
    {
        return D_;
    }

    Real const& D() const
    {
        return D_;
    }

    Species& species()
    {
        return species_;
    }

    Species const& species() const
    {
        return species_;
    }

private:

    Position3 position_;
    Real radius_, D_;
    Species species_;
};

struct ParticleID:
        public Identifier<ParticleID, unsigned long long, int>
{
    typedef Identifier<ParticleID, unsigned long long, int> base_type;

    ParticleID(value_type const& value = value_type(0, 0))
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
    std::size_t operator()(ecell4::ParticleID const& val) const
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

#endif /* __PARTICLE_HPP */
