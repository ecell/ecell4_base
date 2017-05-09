#ifndef ECELL4_COLLISION_HPP
#define ECELL4_COLLISION_HPP

#include "types.hpp"
#include "Real3.hpp"

#include "AABB.hpp"
#include "Cylinder.hpp"
#include "PlanarSurface.hpp"
#include "Sphere.hpp"
#include "Rod.hpp"
#include "Circle.hpp"
#include "Triangle.hpp"
#include "Cone.hpp"
#include "Barycentric.hpp"


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

/* ------ 2D stuff ------ */

//XXX PROPOSAL:
//        define struct Ray and Segment as one-dimensional Shape
//        and use template specialization like
//        - distance_sq<Sphere, Triangle>
//        - test<Sphere, Triangle>
//        - intersect<Ray, Triangle>

Real3 closest_point_point_triangle(const Real3&, const Triangle&);
Real3 closest_point_point_circle(const Real3&, const Circle&);
Real3 closest_point_point_cone(const Real3&, const Cone&);

inline Real distance_sq_point_triangle(const Real3& p, const Triangle& t)
{
    return length_sq(p - closest_point_point_triangle(p, t));
}
inline Real distance_sq_point_circle(const Real3& p, const Circle& c)
{
    return length_sq(p - closest_point_point_circle(p, c));
}
inline Real distance_sq_point_cone(const Real3& p, const Cone& c)
{
    return length_sq(p - closest_point_point_cone(p, c));
}

inline Real distance_sq_sphere_triangle(const Sphere& s, const Triangle& t)
{
    return distance_sq_point_triangle(s.center(), t) - s.radius();
}
inline Real distance_sq_sphere_circle(const Sphere& s, const Circle& c)
{
    return distance_sq_point_circle(s.center(), c) - s.radius();
}
inline Real distance_sq_sphere_cone(const Sphere& s, const Cone& c)
{
    return distance_sq_point_cone(s.center(), c) - s.radius();
}

inline bool test_sphere_triangle(const Sphere& s, const Triangle& t)
{
    return distance_sq_sphere_triangle(s, t) <= s.radius() * s.radius();
}
inline bool test_sphere_circle(const Sphere& s, const Circle& c)
{
    return distance_sq_sphere_circle(s, c) <= s.radius() * s.radius();
}
inline bool test_sphere_cone(const Sphere& s, const Cone& c)
{
    return distance_sq_sphere_cone(s, c) <= s.radius() * s.radius();
}

inline Real distance_sphere_triangle(const Sphere& s, const Triangle& t)
{
    return std::sqrt(distance_sq_sphere_triangle(s, t));
}
inline Real distance_sphere_circle(  const Sphere& s, const Circle& c)
{
    return std::sqrt(distance_sq_sphere_circle(s, c));
}
inline Real distance_sphere_cone(    const Sphere& s, const Cone& c)
{
    return std::sqrt(distance_sq_sphere_cone(s, c));
}

inline Real distance_point_triangle(const Real3& p, const Triangle& t)
{
    return std::sqrt(distance_sq_point_triangle(p, t));
}
inline Real distance_point_circle(  const Real3& p, const Circle& c)
{
    return std::sqrt(distance_sq_point_circle(p, c));
}
inline Real distance_point_cone(    const Real3& p, const Cone& c)
{
    return std::sqrt(distance_sq_point_cone(p, c));
}

bool intersect_segment_triangle(const Real3& p, const Real3& q,
                                const Triangle& t, Barycentric<Real>& b, Real& s);
bool intersect_segment_circle(const Real3& p, const Real3& q,
                              const Circle& c, Real& s);
bool intersect_segment_cone(const Real3& p, const Real3& q,
                            const Cone& c, Real& s);
bool intersect_ray_triangle(const Real3& pos, const Real3& disp,
                            const Triangle& t, Barycentric<Real>& b, Real3& q);
bool intersect_ray_circle(const Real3& pos, const Real3& disp,
                          const Circle& c, Real& t, Real3& q);
bool intersect_ray_cone(const Real3& pos, const Real3& disp,
                        const Cone& c, Real& t, Real3& q);

} // collision

} // ecell4

#endif /* ECELL4_COLLISION_HPP */
