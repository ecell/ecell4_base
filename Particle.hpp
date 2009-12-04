#ifndef PARTICLE_HPP
#define PARTICLE_HPP

#include <ostream>
#if defined(HAVE_TR1_FUNCTIONAL)
#include <tr1/functional>
#elif defined(HAVE_STD_HASH)
#include <functional>
#elif defined(HAVE_BOOST_FUNCTIONAL_HASH_HPP)
#include <boost/functional/hash.hpp>
#endif

#include "Sphere.hpp"

template<typename T_, typename Td_, typename Tsid_>
struct Particle
{
    typedef Sphere<T_> sphere_type;
    typedef Td_ D_type;
    typedef Tsid_ species_id_type;
    typedef typename sphere_type::position_type position_type;
    typedef typename sphere_type::length_type length_type;

    Particle(): sphere_(), species_id_() {}

    Particle(species_id_type const& species_id, sphere_type const& sphere)
        : sphere_(sphere), species_id_(species_id), D_(0.) {}

    position_type& position()
    {
        return sphere_.position();
    }

    position_type const& position() const
    {
        return sphere_.position();
    }

    length_type& radius()
    {
        return sphere_.radius();
    }

    length_type const& radius() const
    {
        return sphere_.radius();
    }

    D_type& D()
    {
        return D_;
    }

    D_type const& D() const
    {
        return D_;
    }

    sphere_type& as_sphere()
    {
        return sphere_;
    }

    sphere_type const& as_sphere() const
    {
        return sphere_;
    }

    species_id_type const& sid() const
    {
        return species_id_;
    }

    species_id_type& sid()
    {
        return species_id_;
    }

    bool operator==(Particle const& rhs) const
    {
        return species_id_ == rhs.sid() && sphere_ == rhs.as_sphere();
    }

    bool operator!=(Particle const& rhs) const
    {
        return !operator==(rhs);
    }

private:
    sphere_type sphere_;
    species_id_type species_id_;
    D_type D_;
};

template<typename Tstrm_, typename Ttraits_, typename T_, typename Td_, typename Tsid_>
inline std::basic_ostream<Tstrm_, Ttraits_>& operator<<(std::basic_ostream<Tstrm_, Ttraits_>& strm, const Particle<T_, Td_, Tsid_>& v)
{
    strm << "Particle(" << v.as_sphere() << ", D=" << v.D() << ", " << v.sid() << ")";
    return strm;
}

#if defined(HAVE_TR1_FUNCTIONAL)
namespace std { namespace tr1 {
#elif defined(HAVE_STD_HASH)
namespace std {
#elif defined(HAVE_BOOST_FUNCTIONAL_HASH_HPP)
namespace boost {
#endif

template<typename T_, typename Td_, typename Tsid_>
struct hash<Particle<T_, Td_, Tsid_> >
{
    typedef Particle<T_, Td_, Tsid_> argument_type;

    std::size_t operator()(argument_type const& val)
    {
        return hash<typename argument_type::position_type>()(val.position()) ^
            hash<typename argument_type::length_type>()(val.radius()) ^
            hash<typename argument_type::D_type>()(val.D()) ^
            hash<typename argument_type::species_id_type>()(val.sid());
    }
};

#if defined(HAVE_TR1_FUNCTIONAL)
} } // namespace std::tr1
#elif defined(HAVE_STD_HASH)
} // namespace std
#elif defined(HAVE_BOOST_FUNCTIONAL_HASH_HPP)
} // namespace boost
#endif

#endif /* PARTICLE_HPP */
