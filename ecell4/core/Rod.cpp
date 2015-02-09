#include "Rod.hpp"
#include "AABB.hpp"
#include "exceptions.hpp"

namespace ecell4
{

Rod::Rod()
    : length_(0.5e-6), radius_(2.0e-6)
{
    ;
}

Rod::Rod(const Real& length, const Real& radius)
    : length_(length), radius_(radius)
{
    //assert(length_ > 0);
    //assert(radius_ > 0)
}

Rod::Rod(const Rod& rhs)
    : length_(rhs.length_), radius_(rhs.radius_)
{
    ;
}

const Real& Rod::lengthX() const
{
    return length_;
}

const Real& Rod::radius() const
{
    return radius_;
}

Real Rod::is_inside(const Real3& pos) const
{
    return distance(pos);
}

Real Rod::distance(const Real3& pos) const
{
    if (pos[0] > length_/2)
        return length(pos - Real3(length_/2, 0, 0)) - radius_;
    if (pos[0] < -length_/2)
        return length(pos - Real3(-length_/2, 0, 0)) - radius_;
    return length(pos - Real3(pos[0], 0, 0)) - radius_;
}

Real3 Rod::draw_position(boost::shared_ptr<RandomNumberGenerator>& rng) const
{
    // The Cylinder Part
    if (rng->uniform(-4*radius_, 3*length_) >= 0)
    {
        const Real x(rng->uniform(-length_/2, length_/2));
        const Real theta(rng->uniform(0, M_PI*2));
        const Real r(sqrtf(rng->uniform(0, radius_*radius_)));
        return Real3(x, r*cos(theta), r*sin(theta));
    }

    // The Terminal Part
    const Real theta(rng->uniform(0, M_PI));
    const Real phi(rng->uniform(0, M_PI));
    const Real r(pow(rng->uniform(0, pow(radius_, 3.0)), 1.0/3.0));
    const Real l(r*sin(phi));
    if (rng->uniform(-1, 1) < 0)
        return Real3(length_/2+l*sin(theta), r*cos(theta), r*cos(phi));
    else
        return Real3(length_/2-l*sin(theta), r*cos(theta), r*cos(phi));
}

RodSurface Rod::surface() const
{
    return RodSurface(length_, radius_);
}

/**
 */

Real clamp(const Real n, const Real min, const Real max)
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


std::pair<bool, Real> intersect_segment_capsule(
    const Real3& p1, const Real3& q1,
    const Real3& p2, const Real3& q2,
    const Real& radius)
{
    Real s, t;
    Real3 c1, c2;
    Real Lsq(closest_point_segment_segment(p1, q1, p2, q2, s, t, c1, c2));
    return std::make_pair(Lsq <= radius * radius, s);
}

/**
 */

bool Rod::test_AABB(const Real3& lower, const Real3& upper) const
{
    const Real3 center(0.0, 0.0, 0.0); //XXX: DEFAULT
    const Real3 axis(1.0, 0.0, 0.0); //XXX: DEFAULT

    const Real3 d(multiply(axis, length_));
    const Real3 p0(center - multiply(axis, length_ * 0.5));
    const Real3 p1(p0 + d);

    const AABB b(lower, upper);
    const AABB e(
        Real3(lower[0] - radius_, lower[1] - radius_, lower[2] - radius_),
        Real3(upper[0] + radius_, upper[1] + radius_, upper[2] + radius_));

    const std::pair<bool, Real> retval(
        e.intersect_ray(p0, d));
    if (!retval.first || retval.second > 1.0)
    {
        return false;
    }

    const Real3 p(p0 + multiply(d, retval.second));
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
        std::pair<bool, Real> retval;
        retval = intersect_segment_capsule(
            p0, p1, b.corner(v), b.corner(v^1), radius_);
        if (retval.first)
        {
            tmin = std::min(retval.second, tmin);
        }
        retval = intersect_segment_capsule(
            p0, p1, b.corner(v), b.corner(v^2), radius_);
        if (retval.first)
        {
            tmin = std::min(retval.second, tmin);
        }
        retval = intersect_segment_capsule(
            p0, p1, b.corner(v), b.corner(v^4), radius_);
        if (retval.first)
        {
            tmin = std::min(retval.second, tmin);
        }

        if (tmin == inf)
        {
            return false;
        }
        return true;
    }

    if ((m & (m - 1)) == 0)
    {
        return true;
    }

    return intersect_segment_capsule(
        p0, p1, b.corner(u^7), b.corner(v), radius_).first;
}

RodSurface::RodSurface()
    : length_(0.5e-6), radius_(2.0e-6)
{
    ;
}

RodSurface::RodSurface(const Real& length, const Real& radius)
    : length_(length), radius_(radius)
{
    //assert(length_ > 0);
    //assert(radius_ > 0)
}

RodSurface::RodSurface(const RodSurface& rhs)
    : length_(rhs.length_), radius_(rhs.radius_)
{
    ;
}

const Real& RodSurface::lengthX() const
{
    return length_;
}

const Real& RodSurface::radius() const
{
    return radius_;
}

Real RodSurface::is_inside(const Real3& pos) const
{
    return distance(pos);
}

Real RodSurface::distance(const Real3& pos) const
{
    if (pos[0] > length_/2)
        return length(pos - Real3(length_/2, 0, 0)) - radius_;
    if (pos[0] < -length_/2)
        return length(pos - Real3(-length_/2, 0, 0)) - radius_;
    return length(pos - Real3(pos[0], 0, 0)) - radius_;
}

Real3 RodSurface::draw_position(boost::shared_ptr<RandomNumberGenerator>& rng) const
{
    // The Cylinder Part
    if (rng->uniform(-2*radius_, length_) >= 0)
    {
        const Real x(rng->uniform(-length_/2, length_/2));
        const Real theta(rng->uniform(0, M_PI*2));
        return Real3(x, radius_*sin(theta), radius_*cos(theta));
    }

    // The Terminal Part
    const Real theta(rng->uniform(0, M_PI));
    const Real phi(rng->uniform(0, M_PI));
    const Real l(radius_*sin(theta));
    if (rng->uniform(-1, 1) < 0)
        return Real3(length_/2+l*sin(theta), l*cos(theta), radius_*cos(phi));
    else
        return Real3(-length_/2-l*sin(theta), l*cos(theta), radius_*cos(phi));
}

Rod RodSurface::inside() const
{
    return Rod(length_, radius_);
}

bool RodSurface::test_AABB(const Real3& lower, const Real3& upper) const
{
    throw NotImplemented("not implemented yet.");
}

} // ecell4
