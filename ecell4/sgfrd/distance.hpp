#ifndef ECELL4_SGFRD_DISTANCE
#define ECELL4_SGFRD_DISTANCE
#include <ecell4/core/collision.hpp>

namespace ecell4
{
namespace sgfrd
{

template<typename shape1, typename shape2>
struct distance_sq
{
    typedef Real result_type;

    Real operator()(const shape1& s1, const shape2& s2) const
    {
        throw ecell4::NotSupported(
                "distance_sq between these shapes is not supported");
    }
};

template<>
struct distance_sq<Real3, ecell4::Triangle>
{
    typedef Real result_type;

    Real operator()(const Real3& p, const ecell4::Triangle& t) const
    {
        return ecell4::collision::distance_sq_point_triangle(p, t);
    }
};

template<>
struct distance_sq<Real3, ecell4::Circle>
{
    typedef Real result_type;

    Real operator()(const Real3& p, const ecell4::Circle& t) const
    {
        return ecell4::collision::distance_sq_point_circle(p, t);
    }
};

template<>
struct distance_sq<Real3, ecell4::Cone>
{
    typedef Real result_type;

    Real operator()(const Real3& p, const ecell4::Cone& t) const
    {
        return ecell4::collision::distance_sq_point_cone(p, t);
    }
};

template<>
struct distance_sq<Real3, ecell4::Cylinder>
{
    typedef Real result_type;

    Real operator()(const Real3& p, const ecell4::Cylinder& t) const
    {
        const Real d = ecell4::collision::distance_point_cylinder(p, t);
        return d*d;
    }
};

template<>
struct distance_sq<Real3, ecell4::Sphere>
{
    typedef Real result_type;

    Real operator()(const Real3& p, const ecell4::Sphere& s) const
    {
        const Real d = length(p - s.center()) - s.radius();
        return d*d;
    }
};

//---------------------------------- distance ----------------------------------

template<typename shape1, typename shape2>
struct distance
{
    typedef Real result_type;

    Real operator()(shape1 s1, shape2 s2) const
    {
        throw ecell4::NotSupported(
                "distance between these shapes is not supported");
    }
};

template<>
struct distance<Real3, ecell4::Triangle>
{
    typedef Real result_type;

    Real operator()(const Real3& p, const ecell4::Triangle& t) const
    {
        return ecell4::collision::distance_point_triangle(p, t);
    }
};

template<>
struct distance<Real3, ecell4::Circle>
{
    typedef Real result_type;

    Real operator()(const Real3& p, const ecell4::Circle& t) const
    {
        return ecell4::collision::distance_point_circle(p, t);
    }
};

template<>
struct distance<Real3, ecell4::Cone>
{
    typedef Real result_type;

    Real operator()(const Real3& p, const ecell4::Cone& t) const
    {
        return ecell4::collision::distance_point_cone(p, t);
    }
};

template<>
struct distance<Real3, ecell4::Cylinder>
{
    typedef Real result_type;

    Real operator()(const Real3& p, const ecell4::Cylinder& t) const
    {
        return ecell4::collision::distance_point_cylinder(p, t);
    }
};

template<>
struct distance<Real3, ecell4::Sphere>
{
    typedef Real result_type;

    Real operator()(const Real3& p, const ecell4::Sphere& s) const
    {
        return length(p - s.center()) - s.radius();
    }
};

} // sgfrd
}// ecell4
#endif// ECELL4_SGFRD_DISTANCE
