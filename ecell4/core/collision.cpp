#include "collision.hpp"


namespace ecell4
{

namespace collision
{

Real distance_sq_point_AABB(const Real3& pos, const AABB& b)
{
    const Real3& upper(b.upper());
    const Real3& lower(b.lower());

    Real Lsq(0.0);
    const unsigned int ndim(3);
    for (unsigned int i(0); i < ndim; ++i)
    {
        const Real& v(pos[i]);
        if (v < lower[i])
        {
            Lsq += pow_2(lower[i] - v);
        }
        else if (v > upper[i])
        {
            Lsq += pow_2(v - upper[i]);
        }
    }
    return Lsq;
}

Real farthest_distance_sq_point_AABB(const Real3& pos, const AABB& b)
{
    const Real3 c(b.center());
    const Real3& lower(b.lower());
    const Real3& upper(b.upper());

    const Real3 q(
        pos[0] > c[0] ? lower[0] : upper[0],
        pos[1] > c[1] ? lower[1] : upper[1],
        pos[2] > c[2] ? lower[2] : upper[2]);
    return length_sq(q - pos);
}

Real distance_point_cylinder(const Real3& pos, const Cylinder& c)
{
    /* First compute the (z,r) components of pos in a coordinate system
     * defined by the vectors unitR and unit_z, where unitR is
     * choosen such that unitR and unit_z define a plane in which
     * pos lies. */
    const Real& half_height(c.half_height());
    const Real& radius(c.radius());
    const std::pair<Real, Real> r_z(c.to_internal(pos));

    /* Then compute distance to cylinder. */
    const Real dz(std::fabs(r_z.second) - half_height);
    const Real dr(r_z.first - radius);

    if (dz > 0)
    {
        // pos is (either) to the right or to the left of the cylinder.
        if (r_z.first > radius)
        {
            // Compute distance to edge.
            return std::sqrt(dz * dz + dr * dr);
        }
        else
        {
            return dz;
        }
    }

    if (dr > radius)
    {
        // pos is somewhere 'parallel' to the cylinder.
        return dr;
    }

    // Inside cylinder.
    return std::max(dr, dz);
}

Real distance_point_capsule(const Real3& pos, const Rod& r)
{
    const Real& half_length(r.half_length());
    const Real& radius(r.radius());
    const Real3 vec(pos - r.origin());

    if (vec[0] > half_length)
    {
        return length(vec - Real3(half_length, 0, 0)) - radius;
    }
    else if (vec[0] < -half_length)
    {
        return length(vec + Real3(half_length, 0, 0)) - radius;
    }
    return length(vec - Real3(vec[0], 0, 0)) - radius;
}

Real closest_point_segment_segment(
    const Real3& p1, const Real3& q1,
    const Real3& p2, const Real3& q2,
    Real& s, Real& t, Real3& c1, Real3& c2)
{
    const Real3 d1(q1 - p1);
    const Real3 d2(q2 - p2);
    const Real3 r(p1 - p2);
    Real a(length_sq(d1));
    Real e(length_sq(d2));
    Real f(dot_product(d2, r));

    if (a <= epsilon && e <= epsilon)
    {
        c1 = p1;
        c2 = p2;
        return length_sq(c1 - c2);
    }

    if (a <= epsilon)
    {
        s = 0.0;
        t = f / e;
        t = clamp(t, 0.0, 1.0);
    }
    else
    {
        Real c = dot_product(d1, r);
        if (e <= epsilon)
        {
            t = 0.0;
            s = clamp(-c / a, 0.0, 1.0);
        }
        else
        {
            Real b = dot_product(d1, d2);
            Real denom = a * e - b * b;
            if (denom != 0.0)
            {
                s = clamp((b * f - c * e)/ denom, 0.0, 1.0);
            }
            else
            {
                s = 0.0;
            }

            t = (b * s + f) / e;
            if (t < 0.0)
            {
                t = 0.0;
                s = clamp(-c / a, 0.0, 1.0);
            }
            else if (t > 1.0)
            {
                t = 1.0;
                s = clamp((b - c) / a, 0.0, 1.0);
            }
        }
    }

    c1 = p1 + d1 * s;
    c2 = p2 + d2 * t;
    return length_sq(c1 - c2);
}

bool test_AABB_AABB(
    const Real3& l1, const Real3& u1,
    const Real3& l2, const Real3& u2)
{
    if (u1[0] < l2[0] || l1[0] > u2[0])
    {
        return false;
    }
    else if (u1[1] < l2[1] || l1[1] > u2[1])
    {
        return false;
    }
    else if (u1[2] < l2[2] || l1[2] > u2[2])
    {
        return false;
    }
    return true;
}

bool test_segment_AABB(
    const Real3& p0, const Real3& p1, const Real3& lower, const Real3& upper)
{
    const Real3 c((upper + lower) * 0.5);
    const Real3 e(upper - c);
    Real3 m(multiply(p1 - p0, 0.5));
    const Real3 d(p1 - m);
    m = m - c;

    Real adx(abs(d[0]));
    if (abs(m[0]) > e[0] + adx)
    {
        return false;
    }
    Real ady(abs(d[1]));
    if (abs(m[1]) > e[1] + ady)
    {
        return false;
    }
    Real adz(abs(d[2]));
    if (abs(m[2]) > e[2] + adz)
    {
        return false;
    }

    adx += epsilon;
    ady += epsilon;
    adz += epsilon;
    if (abs(m[1] * d[2] - m[2] * d[1]) > e[1] * adz + e[2] * ady)
    {
        return false;
    }
    if (abs(m[2] * d[0] - m[0] * d[2]) > e[0] * adz + e[2] * adx)
    {
        return false;
    }
    if (abs(m[0] * d[1] - m[1] * d[0]) > e[0] * ady + e[1] * adx)
    {
        return false;
    }
    return true;
}

bool test_AABB_plane(const AABB& b, const PlanarSurface& p)
{
    const Real3 c(b.center());
    const Real3 e(b.radius());

    const Real3& n(p.normal());
    const Real d(dot_product(p.origin(), n));

    const Real r(e[0] * abs(n[0]) + e[1] * abs(n[1]) + e[2] * abs(n[2]));
    const Real s(dot_product(n, c) - d);
    return (abs(s) <= r);
}

bool test_sphere_AABB(const Sphere& s, const AABB& b)
{
    const Real3& center(s.center());
    const Real& r(s.radius());

    const Real Lsq(distance_sq_point_AABB(center, b));
    return (Lsq <= r * r);
}

bool test_shell_AABB(const SphericalSurface& s, const AABB& b)
{
    const Real r(s.radius());
    const Real rsq(r * r);
    const Real3& center(s.center());

    if (distance_sq_point_AABB(center, b) > rsq)
    {
        return false;
    }
    else if (farthest_distance_sq_point_AABB(center, b) < rsq)
    {
        return false;
    }
    return true;
}

/**
 * intersect_ray_AABB
 * See RTCD p.180;
 */
bool intersect_ray_AABB(
    const Real3& p, const Real3& d, const Real3& lower, const Real3& upper,
    Real& tmin, Real3& q)
{
    tmin = 0.0;
    Real tmax(inf);
    const unsigned int ndim(3);
    for (unsigned int i(0); i < ndim; ++i)
    {
        if (abs(d[i]) < epsilon)
        {
            if (p[i] < lower[i] || p[i] > upper[i])
            {
                return false;
            }
        }
        else
        {
            Real ood = 1.0 / d[i];
            Real t1 = (lower[i] - p[i]) * ood;
            Real t2 = (upper[i] - p[i]) * ood;
            if (t1 > t2)
            {
                const Real tmp(t1);
                t1 = t2;
                t2 = tmp;
            }
            tmin = std::max(tmin, t1);
            tmax = std::min(tmax, t2);
            if (tmin > tmax)
            {
                return false;
            }
        }
    }

    q = p + d * tmin;
    return true;
}

bool intersect_segment_capsule(
    const Real3& p1, const Real3& q1,
    const Real3& p2, const Real3& q2,
    const Real& radius, Real& s)
{
    Real t;
    Real3 c1, c2;
    const Real Lsq(closest_point_segment_segment(p1, q1, p2, q2, s, t, c1, c2));
    return Lsq <= radius * radius;
}

bool intersect_moving_sphere_AABB(
    const Sphere& s, const Real3& d, const AABB& b, Real& t)
{
    const Real3 p0(s.center());
    const Real3 p1(p0 + d);
    const Real& radius(s.radius());
    const Real3& lower(b.lower());
    const Real3& upper(b.upper());

    const AABB e(
        Real3(lower[0] - radius, lower[1] - radius, lower[2] - radius),
        Real3(upper[0] + radius, upper[1] + radius, upper[2] + radius));

    Real3 p;
    if (!intersect_ray_AABB(p0, d, e, t, p) || t > 1.0)
    {
        return false;
    }

    int u(0), v(0);
    if (p[0] < lower[0]) u |= 1;
    if (p[0] > upper[0]) v |= 1;
    if (p[1] < lower[1]) u |= 2;
    if (p[1] > upper[1]) v |= 2;
    if (p[2] < lower[2]) u |= 4;
    if (p[2] > upper[2]) v |= 4;
    const int m(u + v);

    if (m == 7)
    {
        Real tmin(inf);
        if (intersect_segment_capsule(
            p0, p1, b.corner(v), b.corner(v^1), radius, t))
        {
            tmin = std::min(t, tmin);
        }
        if (intersect_segment_capsule(
            p0, p1, b.corner(v), b.corner(v^2), radius, t))
        {
            tmin = std::min(t, tmin);
        }
        if (intersect_segment_capsule(
            p0, p1, b.corner(v), b.corner(v^4), radius, t))
        {
            tmin = std::min(t, tmin);
        }

        if (tmin == inf)
        {
            return false;
        }
        t = tmin;
        return true;
    }

    if ((m & (m - 1)) == 0)
    {
        return true;
    }

    return intersect_segment_capsule(
        p0, p1, b.corner(u^7), b.corner(v), radius, t);
}

} // collision

} // ecell4
