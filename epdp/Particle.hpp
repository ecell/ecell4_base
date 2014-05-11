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

#include <ecell4/core/types.hpp>
#include "Sphere.hpp"
#include "Shape.hpp"

template<typename Tsid_>
struct Particle;

template<typename Tsid_>
inline typename Particle<Tsid_>::shape_type shape(Particle<Tsid_> p)
{
    return p.shape();
}

template<typename Tsid_>
struct Particle
{
    typedef Sphere shape_type;
    typedef ecell4::Real D_type;
    typedef ecell4::Real v_type;	// the drift v has the same type as diffusion constant D for now, may be generalized at a later stage
    typedef Tsid_ species_id_type;
    typedef typename shape_type::position_type position_type;
    typedef typename shape_type::length_type length_type;

    Particle(): shape_(), species_id_(), D_(0.), v_(0.), position_(), radius_(0.)
    {}

    Particle(species_id_type const& species_id, shape_type const& shape,
             D_type const& D)
        : shape_(shape), species_id_(species_id), D_(D), v_(0.), 
            position_(shape.position()), radius_(shape.radius()) {}

    Particle(species_id_type const& species_id, shape_type const& shape,
             D_type const& D, v_type const& v)
        : shape_(shape), species_id_(species_id), D_(D), v_(v),
            position_(shape.position()), radius_(shape.radius()) {}

    position_type& position()
    {
        //return shape_.position();
        return this->position_;
    }

    position_type const& position() const
    {
        //return shape_.position();
        return this->position_;
    }

    length_type& radius()
    {
        //return shape_.radius();
        return this->radius_;
    }

    length_type const& radius() const
    {
        //return shape_.radius();
        return this->radius_;
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

    /*
    shape_type& shape()
    {
        return shape_;
    }
    */

    //shape_type const& shape() const
    shape_type shape() const
    {
        //return shape_;
        return shape_type(this->position(), this->radius());
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
        //return species_id_ == rhs.sid() && shape_ == rhs.shape();
        return (species_id_ == rhs.sid() && radius_ == rhs.radius() &&
                position_ == rhs.position() );
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

    position_type position_;
    length_type radius_;
};

template<typename Tstrm_, typename Ttraits_, typename Tsid_>
inline std::basic_ostream<Tstrm_, Ttraits_>& operator<<(std::basic_ostream<Tstrm_, Ttraits_>& strm, const Particle<Tsid_>& p)
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

template<typename Tsid_>
struct hash<Particle<Tsid_> >
{
    typedef Particle<Tsid_> argument_type;

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
