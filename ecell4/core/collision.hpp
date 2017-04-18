#ifndef ECELL4_COLLISION_HPP
#define ECELL4_COLLISION_HPP

#include "types.hpp"
#include "Real3.hpp"

#include "AABB.hpp"
#include "Cylinder.hpp"
#include "PlanarSurface.hpp"
#include "Sphere.hpp"
#include "Rod.hpp"


namespace ecell4
{

namespace collision
{

inline Real clamp(const Real n, const Real min, const Real max)
{
    if (n < min)
    {
        return min;
    }
    else if (n > max)
    {
        return max;
    }
    return n;
}

Real distance_sq_point_AABB(const Real3& pos, const AABB& b);
Real farthest_distance_sq_point_AABB(const Real3& pos, const AABB& b);
Real distance_point_cylinder(const Real3& pos, const Cylinder& c);
Real distance_point_capsule(const Real3& pos, const Rod& r);

inline Real distance_sq_point_cylinder(const Real3& pos, const Cylinder& c)
{
    const Real L(distance_point_cylinder(pos, c));
    return L * L;
}

Real closest_point_segment_segment(
    const Real3& p1, const Real3& q1,
    const Real3& p2, const Real3& q2,
    Real& s, Real& t, Real3& c1, Real3& c2);

bool test_AABB_AABB(
    const Real3& l1, const Real3& u1,
    const Real3& l2, const Real3& u2);

inline bool test_AABB_AABB(const AABB& b1, const AABB& b2)
{
    return test_AABB_AABB(b1.lower(), b1.upper(), b2.lower(), b2.upper());
}

bool test_segment_AABB(
    const Real3& p0, const Real3& p1, const Real3& lower, const Real3& upper);

inline bool test_segment_AABB(
    const Real3& p0, const Real3& p1, const AABB& b)
{
    return test_segment_AABB(p0, p1, b.lower(), b.upper());
}

bool test_AABB_plane(const AABB& b, const PlanarSurface& p);
bool test_shell_AABB(const SphericalSurface& s, const AABB& b);

inline bool test_shell_AABB(const SphericalSurface& s, const Real3& l, const Real3& u)
{
    return test_shell_AABB(s, AABB(l, u));
}

bool test_sphere_AABB(const Sphere& s, const AABB& b);

inline bool test_sphere_AABB(const Sphere& s, const Real3& l, const Real3& u)
{
    return test_sphere_AABB(s, AABB(l, u));
}

bool intersect_ray_AABB(
    const Real3& p, const Real3& d, const Real3& lower, const Real3& upper,
    Real& tmin, Real3& q);

inline bool intersect_ray_AABB(
    const Real3& p, const Real3& d, const AABB& b,
    Real& tmin, Real3& q)
{
    return intersect_ray_AABB(p, d, b.lower(), b.upper(), tmin, q);
}

bool intersect_segment_capsule(
    const Real3& p1, const Real3& q1,
    const Real3& p2, const Real3& q2,
    const Real& radius, Real& s);

bool intersect_moving_sphere_AABB(
    const Sphere& s, const Real3& d, const AABB& b, Real& t);

} // collision

} // ecell4

#endif /* ECELL4_COLLISION_HPP */
