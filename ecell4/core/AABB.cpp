#include "AABB.hpp"

namespace ecell4
{

Real AABB::distance_sq(const Real3 pos) const
{
    Real Lsq(0.0);
    const unsigned int ndim(3);
    for (unsigned int i(0); i < ndim; ++i)
    {
        const Real& v(pos[i]);
        if (v < lower_[i])
        {
            Lsq += pow_2(lower_[i] - v);
        }
        else if (v > upper_[i])
        {
            Lsq += pow_2(v - upper_[i]);
        }
    }
    return Lsq;
}

Real AABB::distance(const Real3& pos) const
{
    return sqrt(distance_sq(pos));
    // const Real3 closest(
    //     std::min(std::max(pos[0], lower_[0]), upper_[0]),
    //     std::min(std::max(pos[1], lower_[1]), upper_[1]),
    //     std::min(std::max(pos[2], lower_[2]), upper_[2]));
    // const Real distance_plus(length(closest - pos));

    // if (distance_plus > 0)
    // {
    //     return distance_plus;
    // }
    // else
    // {
    //     std::vector<Real> tmp(6);
    //     tmp[0] = lower_[0] - pos[0];
    //     tmp[1] = lower_[1] - pos[1];
    //     tmp[2] = lower_[2] - pos[2];
    //     tmp[3] = pos[0] - upper_[0];
    //     tmp[4] = pos[1] - upper_[1];
    //     tmp[5] = pos[2] - upper_[2];
    //     return *std::max_element(tmp.begin(), tmp.end());
    // }
}

Real3 AABB::draw_position(
    boost::shared_ptr<RandomNumberGenerator>& rng) const
{
    const Real3 pos(
        rng->uniform(lower_[0], upper_[0]),
        rng->uniform(lower_[1], upper_[1]),
        rng->uniform(lower_[2], upper_[2]));
    return pos;
}

Real AABB::farthest_distance_sq(const Real3& pos) const
{
    const Real3 c(center());
    const Real3 q(
        pos[0] > c[0] ? lower_[0] : upper_[0],
        pos[1] > c[1] ? lower_[1] : upper_[1],
        pos[2] > c[2] ? lower_[2] : upper_[2]);
    return distance_sq(q);
}

bool AABB::test_AABB(const Real3& l, const Real3& u) const
{
    if (upper_[0] < l[0] || lower_[0] > u[0])
    {
        return false;
    }
    else if (upper_[1] < l[1] || lower_[1] > u[1])
    {
        return false;
    }
    else if (upper_[2] < l[2] || lower_[2] > u[2])
    {
        return false;
    }
    return true;
}

bool AABB::test_segment(const Real3& p0, const Real3& p1) const
{
    const Real3 c(center());
    const Real3 e(upper_ - c); // radius()
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
    return 1;
}

std::pair<bool, Real> AABB::intersect_ray(const Real3& p, const Real3& d) const
{
    Real tmin(0.0);
    Real tmax(inf);
    const unsigned int ndim(3);
    for (unsigned int i(0); i < ndim; ++i)
    {
        if (abs(d[i]) < epsilon)
        {
            if (p[i] < lower_[i] || p[i] > upper_[i])
            {
                return std::make_pair(false, tmin);
            }
        }
        else
        {
            Real ood = 1.0 / d[i];
            Real t1 = (lower_[i] - p[i]) * ood;
            Real t2 = (upper_[i] - p[i]) * ood;
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
                return std::make_pair(false, tmin);
            }
        }
    }

    // const Real3 q(p + multiply(d, tmin));
    return std::make_pair(true, tmin);
}

}
