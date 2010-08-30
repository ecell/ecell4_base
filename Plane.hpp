#ifndef PLANE_HPP
#define PLANE_HPP

#include <boost/array.hpp>
#include <boost/multi_array.hpp>
#include <utility>
#include <algorithm>
#include "utils/array_helper.hpp"
#include "Shape.hpp"
#include "linear_algebra.hpp"

template<typename T_>
class Plane
{
public:
    typedef T_ value_type;
    typedef Vector3<T_> position_type;
    typedef T_ length_type;

public:
    Plane(position_type const& position = position_type())
        : position_(position),
          units_(array_gen(
            create_vector<position_type>(1., 0., 0.),
            create_vector<position_type>(0., 1., 0.),
            create_vector<position_type>(0., 0., 1.))),
          extent_(array_gen<length_type>(1., 1.)) {}

    template<typename Tarray_>
    Plane(position_type const& position, Tarray_ const& extent)
        : position_(position),
          units_(array_gen(
            create_vector<position_type>(1., 0., 0.),
            create_vector<position_type>(0., 1., 0.),
            create_vector<position_type>(0., 0., 1.)))
    {
        std::copy(boost::begin(extent), boost::end(extent),
                  boost::begin(extent_));
    }

    template<typename Tarray1, typename Tarray2>
    Plane(position_type const& position,
        Tarray1 const& units, Tarray2 const& extent)
        : position_(position)
    {
        std::copy(boost::begin(units), boost::end(units),
                  boost::begin(units_));
        std::copy(boost::begin(extent), boost::end(extent),
                  boost::begin(extent_));
    }

    template<typename Tarray_>
    Plane(position_type const& position,
        position_type const& vx,
        position_type const& vy,
        Tarray_ const& extent = array_gen<length_type>(1., 1.))
        : position_(position), units_(array_gen(vx, vy, cross_product(vx, vy)))
    {
        std::copy(boost::begin(extent), boost::end(extent),
                  boost::begin(extent_));
    }

    Plane(position_type const& position,
        position_type const& vx,
        position_type const& vy,
        length_type const& lx,
        length_type const& ly)
        : position_(position), units_(array_gen(vx, vy, cross_product(vx, vy))),
          extent_(array_gen<length_type>(lx, ly)) {}

    position_type const& position() const
    {
        return position_;
    }

    position_type& position()
    {
        return position_;
    }

    position_type const& unit_x() const
    {
        return units_[0];
    }

    position_type& unit_x()
    {
        return units_[0];
    }

    position_type const& unit_y() const
    {
        return units_[1];
    }

    position_type& unit_y()
    {
        return units_[1];
    }

    position_type const& unit_z() const
    {
        return units_[2];
    }

    position_type& unit_z()
    {
        return units_[2];
    }

    boost::array<position_type, 3> const& units() const
    {
        return units_;
    }

    boost::array<position_type, 3>& units()
    {
        return units_;
    }

    length_type const& Lx() const
    { 
        return extent_[0];
    }

    length_type& Lx()
    {
        return extent_[0];
    }

    length_type const& Ly() const
    {
        return extent_[1];
    }

    length_type& Ly()
    {
        return extent_[1];
    }

    boost::array<length_type, 2> const& extent() const
    {
        return extent_;
    }

    boost::array<length_type, 2>& extent()
    {
        return extent_;
    }

    bool operator==(const Plane& rhs) const
    {
        return position_ == rhs.position_ && units_ == rhs.units_ &&
               extent_ == rhs.extent_;
    }

    bool operator!=(const Plane& rhs) const
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

protected:
    position_type position_;
    boost::array<position_type, 3> units_;
    boost::array<length_type, 2> extent_;
};

template<typename T_>
inline boost::array<typename Plane<T_>::length_type, 3>
to_internal(Plane<T_> const& obj, typename Plane<T_>::position_type const& pos)
{
    typedef typename Plane<T_>::position_type position_type;
    position_type pos_vector(subtract(pos, obj.position()));

    return array_gen<typename Plane<T_>::length_type>(
        dot_product(pos_vector, obj.unit_x()),
        dot_product(pos_vector, obj.unit_y()),
        dot_product(pos_vector, obj.unit_z()));
}

template<typename T_>
inline std::pair<typename Plane<T_>::position_type,
                 typename Plane<T_>::length_type>
projected_point(Plane<T_> const& obj, typename Plane<T_>::position_type const& pos)
{
    boost::array<typename Plane<T_>::length_type, 3> x_y_z(to_internal(obj, pos));
    return std::make_pair(
        add(add(obj.position(), multiply(obj.unit_x(), x_y_z[0])),
            multiply(obj.unit_y(), x_y_z[1])),
        x_y_z[2]);
}

template<typename T_>
inline typename Plane<T_>::length_type
distance(Plane<T_> const& obj, typename Plane<T_>::position_type const& pos)
{
    typedef typename Plane<T_>::length_type length_type;
    boost::array<length_type, 3> const x_y_z(to_internal(obj, pos));

    length_type const dx(subtract(abs(x_y_z[0]), obj.extent()[0]));
    length_type const dy(subtract(abs(x_y_z[1]), obj.extent()[1]));

    if (dx < 0 && dy < 0) {
        // Projected point of pos is on the plane.
        // Probably an infinite plane anyway.
        return x_y_z[2];
    }

    if (dx > 0)
    {
        if (dy > 0)
        {
            // Far away from plane.
            return std::sqrt(gsl_pow_2(dx) + gsl_pow_2(dy) +
                             gsl_pow_2(x_y_z[2]));
        }
        else
        {
            return std::sqrt(gsl_pow_2(dx) + gsl_pow_2(x_y_z[2]));
        }
    }
    else
    {
        if (dy > 0)
        {
            return std::sqrt(gsl_pow_2(dy) + gsl_pow_2(x_y_z[2]));
        }
        else
        {
            // Already tested above.
            return x_y_z[2];
        }
    }
}

template<typename T, typename Trng>
inline typename Plane<T>::position_type
random_position(Plane<T> const& shape, Trng& rng)
{
    typedef typename Plane<T>::length_type length_type;

    // -1 < rng() < 1. See for example PlanarSurface.hpp.
    return add(
        shape.position(),
        add(multiply(shape.units()[0], shape.extent()[0] * rng()),
            multiply(shape.units()[1], shape.extent()[1] * rng())));
}

template<typename T>
inline Plane<T> const& shape(Plane<T> const& shape)
{
    return shape;
}

template<typename T>
inline Plane<T>& shape(Plane<T>& shape)
{
    return shape;
}

template<typename T_>
struct is_shape<Plane<T_> >: public boost::mpl::true_ {};

template<typename T_>
struct shape_position_type<Plane<T_> >
{
    typedef typename Plane<T_>::position_type type;
};

template<typename Tstrm_, typename Ttraits_, typename T_>
inline std::basic_ostream<Tstrm_, Ttraits_>& operator<<(std::basic_ostream<Tstrm_, Ttraits_>& strm,
        const Plane<T_>& v)
{
    strm << "{" << v.position() <<  ", " << v.unit_x() << ", " << v.unit_y() << "," << v.Lx() << ", " << v.Ly() << "}";
    return strm;
}


#if defined(HAVE_TR1_FUNCTIONAL)
namespace std { namespace tr1 {
#elif defined(HAVE_STD_HASH)
namespace std {
#elif defined(HAVE_BOOST_FUNCTIONAL_HASH_HPP)
namespace boost {
#endif

template<typename T_>
struct hash<Plane<T_> >
{
    typedef Plane<T_> argument_type;

    std::size_t operator()(argument_type const& val)
    {
        return hash<typename argument_type::position_type>()(val.position()) ^
            hash<typename argument_type::position_type>()(val.unit_x()) ^
            hash<typename argument_type::position_type>()(val.unit_y()) ^
            hash<typename argument_type::length_type>()(val.extent()[0]) ^
            hash<typename argument_type::length_type>()(val.extent()[1]);
    }
};

#if defined(HAVE_TR1_FUNCTIONAL)
} } // namespace std::tr1
#elif defined(HAVE_STD_HASH)
} // namespace std
#elif defined(HAVE_BOOST_FUNCTIONAL_HASH_HPP)
} // namespace boost
#endif

#endif /* PLANE_HPP */
