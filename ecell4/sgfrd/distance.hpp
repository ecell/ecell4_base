#ifndef ECELL4_SGFRD_DISTANCE
#define ECELL4_SGFRD_DISTANCE
#include <ecell4/core/Segment.hpp>
#include <ecell4/core/collision.hpp>
#include <boost/type_traits.hpp>
#include <boost/mpl/and.hpp>
#include <boost/mpl/or.hpp>
#include <boost/format.hpp>

namespace ecell4
{
namespace sgfrd
{

template<typename shape1, typename shape2>
struct distance_sq_impl
{
    typedef Real result_type;

    Real operator()(const shape1& s1, const shape2& s2) const
    {
        throw ecell4::NotSupported(
                "distance_sq for this specific shapes is not supported");
    }
};

template<>
struct distance_sq_impl<Real3, ecell4::Triangle>
{
    typedef Real result_type;

    Real operator()(const Real3& p, const ecell4::Triangle& t) const
    {
        return ecell4::collision::distance_sq_point_triangle(p, t);
    }
};

template<>
struct distance_sq_impl<Real3, ecell4::Circle>
{
    typedef Real result_type;

    Real operator()(const Real3& p, const ecell4::Circle& t) const
    {
        return ecell4::collision::distance_sq_point_circle(p, t);
    }
};

template<>
struct distance_sq_impl<Real3, ecell4::Cone>
{
    typedef Real result_type;

    Real operator()(const Real3& p, const ecell4::Cone& t) const
    {
        return ecell4::collision::distance_sq_point_cone(p, t);
    }
};

template<>
struct distance_sq_impl<Real3, ecell4::Cylinder>
{
    typedef Real result_type;

    Real operator()(const Real3& p, const ecell4::Cylinder& t) const
    {
        const Real d = ecell4::collision::distance_point_cylinder(p, t);
        return d*d;
    }
};

template<>
struct distance_sq_impl<Real3, ecell4::Sphere>
{
    typedef Real result_type;

    Real operator()(const Real3& p, const ecell4::Sphere& s) const
    {
        const Real d = length(p - s.center()) - s.radius();
        return d*d;
    }
};

template<>
struct distance_sq_impl<Real3, ecell4::Segment>
{
    typedef Real result_type;

    Real operator()(const Real3& p, const ecell4::Segment& s) const
    {
        const Real3 ab = s.stop() - s.start();
        const Real3 ac = p - s.start();
        const Real3 bc = p - s.stop();
        const Real dot = dot_product(ac, ab);
        if(dot <= 0.0) return length_sq(ac);
        const Real len = length_sq(ab);
        if(dot >= len) return length_sq(bc);
        return length_sq(ac) - (dot * dot) / len;
    }
};

//---------------------------------- distance ----------------------------------

template<typename shape1, typename shape2>
struct distance_impl
{
    typedef Real result_type;

    Real operator()(shape1 s1, shape2 s2) const
    {
        throw ecell4::NotSupported(
                "distance for this specific shapes is not supported");
    }
};

template<>
struct distance_impl<Real3, ecell4::Triangle>
{
    typedef Real result_type;

    Real operator()(const Real3& p, const ecell4::Triangle& t) const
    {
        return ecell4::collision::distance_point_triangle(p, t);
    }
};

template<>
struct distance_impl<Real3, ecell4::Circle>
{
    typedef Real result_type;

    Real operator()(const Real3& p, const ecell4::Circle& t) const
    {
        return ecell4::collision::distance_point_circle(p, t);
    }
};

template<>
struct distance_impl<Real3, ecell4::Cone>
{
    typedef Real result_type;

    Real operator()(const Real3& p, const ecell4::Cone& t) const
    {
        return ecell4::collision::distance_point_cone(p, t);
    }
};

template<>
struct distance_impl<Real3, ecell4::Cylinder>
{
    typedef Real result_type;

    Real operator()(const Real3& p, const ecell4::Cylinder& t) const
    {
        return ecell4::collision::distance_point_cylinder(p, t);
    }
};

template<>
struct distance_impl<Real3, ecell4::Sphere>
{
    typedef Real result_type;

    Real operator()(const Real3& p, const ecell4::Sphere& s) const
    {
        return length(p - s.center()) - s.radius();
    }
};

template<>
struct distance_impl<Real3, ecell4::Segment>
{
    typedef Real result_type;

    Real operator()(const Real3& p, const ecell4::Segment& s) const
    {
        return std::sqrt(distance_sq_impl<Real3, ecell4::Segment>()(p, s));
    }
};

// ----------------------------------------------------------------------------

template<typename T1, typename T2>
inline typename boost::enable_if<boost::mpl::and_<
    boost::is_base_of<ecell4::Shape, T1>, boost::is_base_of<ecell4::Shape, T2>
    >, Real>::type
distance(const T1& shape1, const T2& shape2)
{
    return distance_impl<T1, T2>()(shape1, shape2);
}

template<typename T>
inline typename boost::enable_if<boost::is_base_of<ecell4::Shape, T>, Real>::type
distance(const T& shape, const Real3& pos)
{
    return distance_impl<Real3, T>()(pos, shape);
}

template<typename T>
inline typename boost::enable_if<boost::is_base_of<ecell4::Shape, T>, Real>::type
distance(const Real3& pos, const T& shape)
{
    return distance_impl<Real3, T>()(pos, shape);
}

inline Real distance(const Real3& lhs, const Real3& rhs)
{
    return length(lhs - rhs);
}

template<typename T1, typename T2>
inline typename boost::enable_if<boost::mpl::and_<
    boost::is_base_of<ecell4::Shape, T1>, boost::is_base_of<ecell4::Shape, T2>
    >, Real>::type
distance_sq(const T1& shape1, const T2& shape2)
{
    return distance_sq_impl<T1, T2>()(shape1, shape2);
}

template<typename T>
inline typename boost::enable_if<boost::is_base_of<ecell4::Shape, T>, Real>::type
distance_sq(const T& shape, const Real3& pos)
{
    return distance_sq_impl<Real3, T>()(pos, shape);
}

template<typename T>
inline typename boost::enable_if<boost::is_base_of<ecell4::Shape, T>, Real>::type
distance_sq(const Real3& pos, const T& shape)
{
    return distance_sq_impl<Real3, T>()(pos, shape);
}

inline Real distance_sq(const Real3& lhs, const Real3& rhs)
{
    return length_sq(lhs - rhs);
}

} // sgfrd
}// ecell4
#endif// ECELL4_SGFRD_DISTANCE
