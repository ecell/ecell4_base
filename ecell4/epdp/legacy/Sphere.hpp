#ifndef SPHERE_HPP
#define SPHERE_HPP

#include <ostream>
#include "Vector3.hpp"
#include "Position3Type.hpp"
#include "Shape.hpp"

class Sphere;
template<typename Tstrm_, typename Ttraits_>
inline std::basic_ostream<Tstrm_, Ttraits_>& operator<<(std::basic_ostream<Tstrm_, Ttraits_>& strm,
        const Sphere& v);

class Sphere
{
public:
    /*
    typedef T_ value_type;
    typedef Vector3<T_> position_type;
    typedef T_ length_type;
    */
    typedef ecell4::Position3 position_type;
    typedef position_type::value_type value_type;
    typedef position_type::value_type length_type;

public:
    Sphere()
        : position_(), radius_(0) {}

    Sphere(const position_type& position, const length_type& radius)
        : position_(position), radius_(radius) {}

    bool operator==(const Sphere& rhs) const
    {
        return position_ == rhs.position_ && radius_ == rhs.radius_;
    }

    bool operator!=(const Sphere& rhs) const
    {
        return !operator==(rhs);
    }

    position_type const& position() const
    {
        return position_;
    }

    position_type& position()
    {
        return position_;
    }

    length_type const& radius() const
    {
        return radius_;
    }

    length_type& radius()
    {
        return radius_;
    }

    std::string show(int precision)
    {
        std::ostringstream strm;
        strm.precision(precision);
        strm << *this;
        return strm.str();
    }

private:
    position_type position_;
    length_type radius_;
};

template<typename Tstrm_, typename Ttraits_>
inline std::basic_ostream<Tstrm_, Ttraits_>& operator<<(std::basic_ostream<Tstrm_, Ttraits_>& strm,
        const Sphere& v)
{
    strm << "{" << v.position() <<  ", " << v.radius() << "}";
    return strm;
}


inline std::pair<Sphere::position_type,
                 Sphere::length_type>
projected_point(Sphere const& obj,
                Sphere::position_type const& pos)
{
    // Todo. If we ever need it.
    // The projection of a point on a sphere.
    return std::make_pair(Sphere::position_type(),
                          Sphere::length_type());
}


inline Sphere::length_type
distance(Sphere const& obj, Sphere::position_type const& pos)
{
    return distance(pos, obj.position()) - obj.radius();
}


inline Sphere const& shape(Sphere const& shape)
{
    return shape;
}


inline Sphere& shape(Sphere& shape)
{
    return shape;
}

template<typename Trng>
inline Sphere::position_type
random_position(Sphere const& shape, Trng& rng)
{
    return add(shape.position(),
                create_vector<Sphere::position_type>(
                    shape.radius() * rng(),
                    shape.radius() * rng(),
                    shape.radius() * rng())); 
}

template<>
struct is_shape<Sphere>: public boost::mpl::true_ {};

template<>
struct shape_position_type<Sphere> {
    typedef Sphere::position_type type;
};

template<>
struct shape_position_type<const Sphere> {
    typedef Sphere::position_type type;
};

template<>
struct shape_length_type<Sphere> {
    typedef Sphere::length_type type;
};


inline typename shape_length_type<Sphere>::type const& shape_size(Sphere const& shape)
{
    return shape.radius();
} 


inline typename shape_length_type<Sphere>::type& shape_size(Sphere &shape)
{
    return shape.radius();
} 

#if defined(HAVE_TR1_FUNCTIONAL)
namespace std { namespace tr1 {
#elif defined(HAVE_STD_HASH)
namespace std {
#elif defined(HAVE_BOOST_FUNCTIONAL_HASH_HPP)
namespace boost {
#endif

template<>
struct hash<Sphere>
{
    typedef Sphere argument_type;

    std::size_t operator()(argument_type const& val)
    {
        return hash<typename argument_type::position_type>()(val.position()) ^
            hash<typename argument_type::length_type>()(val.radius());
    }
};

#if defined(HAVE_TR1_FUNCTIONAL)
} } // namespace std::tr1
#elif defined(HAVE_STD_HASH)
} // namespace std
#elif defined(HAVE_BOOST_FUNCTIONAL_HASH_HPP)
} // namespace boost
#endif

#endif /* SPHERE_HPP */
