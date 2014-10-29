#ifndef CYLINDER_HPP
#define CYLINDER_HPP

#include <ostream>
#include <cmath>
#include "Vector3.hpp"
#include "Position3Type.hpp"
#include "Shape.hpp"

class Cylinder;

template<typename Tstrm_>
inline std::basic_ostream<Tstrm_>& operator<<(std::basic_ostream<Tstrm_>& strm,
        const Cylinder& v);

// Todo. Make sure cylinder is never larger than 1 cellsize or something.  

class Cylinder
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
    Cylinder()
        : position_(), radius_(0), unit_z_(), half_length_(0) {}

    Cylinder(position_type const& position, length_type const& radius,
             position_type const& unit_z, length_type const& half_length )
        : position_(position), radius_(radius), unit_z_(unit_z),
          half_length_(half_length) {}

    bool operator==(const Cylinder& rhs) const
    {
        return position_ == rhs.position() && radius_ == rhs.radius() && unit_z_ == rhs.unit_z() && half_length_ == rhs.half_length();
    }

    bool operator!=(const Cylinder& rhs) const
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

    position_type const& unit_z() const
    {
        return unit_z_;
    }

    position_type& unit_z()
    {
        return unit_z_;
    }

    length_type const& half_length() const
    {
        return half_length_;
    }

    length_type& half_length()
    {
        return half_length_;
    }

    std::string show(int precision)
    {
        std::ostringstream strm;
        strm.precision(precision);
        strm << *this;
        return strm.str();
    }

private:
    position_type position_; // centre.
    length_type radius_;
    position_type unit_z_; // Z-unit_z. should be normalized.
    length_type half_length_;
};

template<typename Tstrm_>
inline std::basic_ostream<Tstrm_>& operator<<(std::basic_ostream<Tstrm_>& strm,
        const Cylinder& v)
{
    strm << "{" << v.position() <<  ", " << v.radius() << ", " << v.unit_z() << ", " << v.half_length() << "}";
    return strm;
}


inline std::pair<Cylinder::length_type,
                 Cylinder::length_type>
to_internal(Cylinder const& obj, Cylinder::position_type const& pos)
{
    // Return pos relative to position of cylinder. 
    typedef Cylinder::position_type position_type;
    typedef Cylinder::length_type length_type;

    const position_type pos_vector(subtract(pos, obj.position()));
    // z can be < 0
    const length_type z(dot_product(pos_vector, obj.unit_z()));
    // r is always >= 0
    const length_type r(length(pos_vector - multiply(obj.unit_z(), z)));

    return std::make_pair(r, z);
}


inline std::pair<Cylinder::position_type,
                 Cylinder::length_type>
projected_point(Cylinder const& obj,
                Cylinder::position_type const& pos)
{
    typedef Cylinder::length_type length_type;

    // The projection lies on the z-axis.
    std::pair<length_type, length_type> r_z(to_internal(obj, pos));
    return std::make_pair(
        add(obj.position(), multiply(obj.unit_z(), r_z.second)),
        r_z.first);
}


inline Cylinder::length_type
distance(Cylinder const& obj,
                Cylinder::position_type const& pos)
{
    typedef Cylinder::position_type position_type;
    typedef Cylinder::length_type length_type;

    /* First compute the (z,r) components of pos in a coordinate system 
     * defined by the vectors unitR and unit_z, where unitR is
     * choosen such that unitR and unit_z define a plane in which
     * pos lies. */
    const std::pair<length_type, length_type> r_z(to_internal(obj, pos));

    /* Then compute distance to cylinder. */
    const length_type dz(std::fabs(r_z.second) - obj.half_length());
    const length_type dr(r_z.first - obj.radius());
    length_type distance;
    if (dz > 0)
    {
        // pos is (either) to the right or to the left of the cylinder.
        if (r_z.first > obj.radius())
        {
            // Compute distance to edge.
            distance = std::sqrt( dz * dz + dr * dr );
        }
        else
        {
            distance = dz;
        }
    }
    else
    {
        if (dr > obj.radius())
        {
            // pos is somewhere 'parallel' to the cylinder.
            distance = dr;
        }
        else
        {
            // Inside cylinder. 
            distance = std::max(dr, dz);
        }
    }
    return distance;
}

template<typename Trng>
inline Cylinder::position_type
random_position(Cylinder const& shape, Trng& rng)
{
    // -1 < rng() < 1. See for example CylindricalSurface.hpp.
    return add(shape.position(),
               multiply(shape.unit_z(), rng() * shape.half_length()));
}


inline Cylinder const& shape(Cylinder const& shape)
{
    return shape;
}


inline Cylinder& shape(Cylinder& shape)
{
    return shape;
}

template<>
struct is_shape<Cylinder>: public boost::mpl::true_ {};

template<>
struct shape_position_type<Cylinder>
{
    typedef typename Cylinder::position_type type;
};

template<>
struct shape_position_type<const Cylinder>
{
    typedef typename Cylinder::position_type type;
};

template<>
struct shape_length_type<Cylinder> {
    typedef typename Cylinder::length_type type;
};


inline typename shape_length_type<Cylinder>::type const& shape_size(Cylinder const& shape)
{
    return shape.radius();
} 


inline typename shape_length_type<Cylinder>::type& shape_size(Cylinder& shape)
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
struct hash<Cylinder>
{
    typedef Cylinder argument_type;

    std::size_t operator()(argument_type const& val)
    {
        return hash<typename argument_type::position_type>()(val.position()) ^
            hash<typename argument_type::length_type>()(val.radius()) ^
            hash<typename argument_type::position_type>()(val.unit_z()) ^
            hash<typename argument_type::length_type>()(val.half_length());
    }
};

#if defined(HAVE_TR1_FUNCTIONAL)
} } // namespace std::tr1
#elif defined(HAVE_STD_HASH)
} // namespace std
#elif defined(HAVE_BOOST_FUNCTIONAL_HASH_HPP)
} // namespace boost
#endif

#endif /* CYLINDER_HPP */
