#include "AABB.hpp"
#include "collision.hpp"


namespace ecell4
{

Real AABB::distance_sq(const Real3 pos) const
{
    return collision::distance_sq_point_AABB(pos, *this);
}

Real AABB::distance(const Real3& pos) const
{
    return sqrt(distance_sq(pos));
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

bool AABB::test_AABB(const Real3& l, const Real3& u) const
{
    return collision::test_AABB_AABB(lower_, upper_, l, u);
}

bool AABB::test_segment(const Real3& p0, const Real3& p1) const
{
    return collision::test_segment_AABB(p0, p1, lower_, upper_);
}

std::pair<bool, Real> AABB::intersect_ray(const Real3& p, const Real3& d) const
{
    Real tmin;
    Real3 q;
    const bool retval(collision::intersect_ray_AABB(p, d, lower_, upper_, tmin, q));
    return std::make_pair(retval, tmin);
}

}
