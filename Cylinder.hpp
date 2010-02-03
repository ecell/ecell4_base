
#ifndef CYLINDER_HPP
#define CYLINDER_HPP

#include <ostream>
#include <cmath>
#include "Vector3.hpp"
#include "Shape.hpp"

// Todo. Make sure cylinder is never larger than 1 cellsize or something.  
template<typename T_>
class Cylinder
{
public:
    typedef T_ value_type;
    typedef Vector3<T_> position_type;
    typedef T_ length_type;

public:
    Cylinder()
        : position_(), radius_(0), ori_(), size_(0) {}

    Cylinder(position_type const& position, length_type const& radius,
             position_type const& ori_z, length_type const& size )
        : position_(position), radius_(radius), ori_(ori_z),
          size_(size) {}

    bool operator==(const Cylinder& rhs) const
    {
        return position_ == rhs.position() && radius_ == rhs.radius() && ori_ == rhs.orientation() && size_ == rhs.size();
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

    position_type const& orientation() const
    {
        return ori_;
    }

    position_type& orientation()
    {
        return ori_;
    }

    length_type const& size() const
    {
        return size_;
    }

    length_type& size()
    {
        return size_;
    }

private:
    position_type position_; // centre.
    length_type radius_;
    position_type ori_; // Z-orientation. should be normalized.
    length_type size_; // half length.
};

template<typename Tstrm_, typename T_>
inline std::basic_ostream<Tstrm_>& operator<<(std::basic_ostream<Tstrm_>& strm,
        const Cylinder<T_>& v)
{
    strm << "{" << v.position() <<  ", " << v.radius() << ", " << v.orientation() << ", " << v.size() << "}";
    return strm;
}

template<typename T_>
inline typename Cylinder<T_>::length_type
distance(Cylinder<T_> const& obj,
         typename Cylinder<T_>::position_type const& pos)
{
    typedef typename Cylinder<T_>::position_type position_type;
    typedef typename Cylinder<T_>::length_type length_type;

    /* First compute the (z,r) components of pos in a coordinate system 
     * defined by the vectors unitR and orientation, where unitR is
     * choosen such that unitR and orientation define a plane in which
     * pos lies. */
    const position_type v = pos - obj.position();
    const length_type z = dot_product(v, obj.orientation()); // Can be <0.
    const position_type pos_z = multiply(obj.orientation(), z);
    const position_type pos_r = v - pos_z;
    const length_type r = length(pos_r);

    /* Then compute distance to cylinder. */
    const length_type dz = std::fabs(z) - obj.size();
    const length_type dr = r - obj.radius();
    length_type distance;
    if (dz > 0)
    {
        // pos is (either) to the right or to the left of the cylinder.
        if (r > obj.radius())
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
            // pos is somewhere 'parellel' to the cylinder.
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

template<typename T_>
inline typename Cylinder<T_>::length_type
distance_cyclic(Cylinder<T_> const& obj, typename Cylinder<T_>::position_type const& pos, typename Cylinder<T_>::length_type const& world_size)
{
    return distance(obj, pos);
}

template<typename T_>
inline Cylinder<T_> const& shape(Cylinder<T_> const& shape)
{
    return shape;
}

template<typename T_>
inline Cylinder<T_>& shape(Cylinder<T_>& shape)
{
    return shape;
}

#if defined(HAVE_TR1_FUNCTIONAL)
namespace std { namespace tr1 {
#elif defined(HAVE_STD_HASH)
namespace std {
#elif defined(HAVE_BOOST_FUNCTIONAL_HASH_HPP)
namespace boost {
#endif

template<typename T_>
struct hash<Cylinder<T_> >
{
    typedef Cylinder<T_> argument_type;

    std::size_t operator()(argument_type const& val)
    {
        return hash<typename argument_type::position_type>()(val.position()) ^
            hash<typename argument_type::length_type>()(val.radius()) ^
            hash<typename argument_type::position_type>()(val.orientation()) ^
            hash<typename argument_type::length_type>()(val.size());
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
