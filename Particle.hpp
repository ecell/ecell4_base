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
#include "Shape.hpp"

template<typename T_, typename Td_, typename Tsid_>
struct Particle
{
    typedef Sphere<T_> shape_type;
    typedef Td_ D_type;
    typedef Td_ v_type;	// the drift v has the same type as diffusion constant D for now, may be generalized at a later stage
    typedef Tsid_ species_id_type;
    typedef typename shape_type::position_type position_type;
    typedef typename shape_type::length_type length_type;

    Particle(): shape_(), species_id_(), D_(0.), v_(0.) {}

    Particle(species_id_type const& species_id, shape_type const& shape,
             D_type const& D)
        : shape_(shape), species_id_(species_id), D_(D), v_(0.){}

    Particle(species_id_type const& species_id, shape_type const& shape,
             D_type const& D, v_type const& v)
        : shape_(shape), species_id_(species_id), D_(D), v_(v) {}

    position_type& position()
    {
        return shape_.position();
    }

    position_type const& position() const
    {
        return shape_.position();
    }

    length_type& radius()
    {
        return shape_.radius();
    }

    length_type const& radius() const
    {
        return shape_.radius();
    }

    D_type& D()
    {
        return D_;
    }

    D_type const& D() const
    {
        return D_;
    }
    
    v_type& v()
    {
        return v_;
    }

    v_type const& v() const
    {
        return v_;
    }

    shape_type& shape()
    {
        return shape_;
    }

    shape_type const& shape() const
    {
        return shape_;
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
        return species_id_ == rhs.sid() && shape_ == rhs.shape();
    }

    bool operator!=(Particle const& rhs) const
    {
        return !operator==(rhs);
    }

    std::string show(int precision)
    {
        std::ostringstream strm;
        strm.precision(precision);
        strm << *this;
        return strm.str();
    }

private:
    shape_type shape_;
    species_id_type species_id_;
    D_type D_;
    v_type v_;
};

template<typename Tstrm_, typename Ttraits_, typename T_, typename Td_, typename Tsid_>
inline std::basic_ostream<Tstrm_, Ttraits_>& operator<<(std::basic_ostream<Tstrm_, Ttraits_>& strm, const Particle<T_, Td_, Tsid_>& p)
{
    strm << "Particle(" << p.shape() << ", D=" << p.D() << ", v=" << p.v() << ", " << p.sid() << ")";
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
