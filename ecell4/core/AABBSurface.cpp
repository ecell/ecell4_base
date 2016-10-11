#include "AABB.hpp"
#include "AABBSurface.hpp"
#include "collision.hpp"


namespace ecell4
{

Real AABBSurface::distance_sq(const Real3 pos) const
{
    if(this->_is_inside(pos))
    {
        const Real3 dupper = upper_ - pos;
        const Real3 dlower = lower_ - pos;
        const Real du = std::min(std::min(dupper[0], dupper[1]), dupper[2]);
        const Real dl = std::min(std::min(dlower[0], dlower[1]), dlower[2]);
        const Real dmin = std::min(du, dl);
        return dmin * dmin;
    }
    else
    {
        return collision::distance_sq_point_AABB(pos, AABB(lower_, upper_));
    }
}

Real AABBSurface::distance(const Real3& pos) const
{
    if(this->_is_inside(pos))
    {
        const Real3 dupper = upper_ - pos;
        const Real3 dlower = lower_ - pos;
        const Real du = std::min(std::min(dupper[0], dupper[1]), dupper[2]);
        const Real dl = std::min(std::min(dlower[0], dlower[1]), dlower[2]);
        return std::min(du, dl);
    }
    else
    {
        return sqrt(collision::distance_sq_point_AABB(pos, AABB(lower_, upper_)));
    }
}

Real3 AABBSurface::draw_position(
    boost::shared_ptr<RandomNumberGenerator>& rng) const
{
    const Real Sxy = (upper_[0] - lower_[0]) * (upper_[1] - lower_[1]);
    const Real Syz = (upper_[1] - lower_[1]) * (upper_[2] - lower_[2]);
    const Real Szx = (upper_[1] - lower_[1]) * (upper_[0] - lower_[0]);
    const Real Stot = (Sxy + Syz + Szx) * 2;
    Real rnd = rng->uniform(0., Stot);

    if((rnd -= Sxy) < 0.)
    {
        return Real3(rng->uniform(lower_[0], upper_[0]),
                     rng->uniform(lower_[1], upper_[1]),
                     lower_[2]);
    }
    else if((rnd -= Sxy) < 0.)
    {
        return Real3(rng->uniform(lower_[0], upper_[0]),
                     rng->uniform(lower_[1], upper_[1]),
                     upper_[2]);
    }
    else if((rnd -= Syz) < 0.)
    {
        return Real3(lower_[0],
                     rng->uniform(lower_[1], upper_[1]),
                     rng->uniform(lower_[2], upper_[2]));
    }
    else if((rnd -= Syz) < 0.)
    {
        return Real3(upper_[0],
                     rng->uniform(lower_[1], upper_[1]),
                     rng->uniform(lower_[2], upper_[2]));
    }
    else if((rnd -= Szx) < 0.)
    {
        return Real3(rng->uniform(lower_[0], upper_[0]),
                     lower_[1],
                     rng->uniform(lower_[2], upper_[2]));
    }
    else if((rnd -= Szx) < 0.)
    {
        return Real3(rng->uniform(lower_[0], upper_[0]),
                     upper_[1],
                     rng->uniform(lower_[2], upper_[2]));
    }
    else
    {
        throw std::logic_error("invalid random number");
    }
}

bool AABBSurface::test_AABB(const Real3& l, const Real3& u) const
{// same as AABB case?
    return collision::test_AABB_AABB(lower_, upper_, l, u);
}

bool AABBSurface::test_segment(const Real3& p0, const Real3& p1) const
{
    if(this->_is_inside(p0) && this->_is_inside(p1))
    {
        return false;
    }
    else if(this->_is_inside(p0) || this->_is_inside(p1))
    {
        return true;
    }
    else
    {
        return collision::test_segment_AABB(p0, p1, lower_, upper_);
    }
}

std::pair<bool, Real> AABBSurface::intersect_ray(const Real3& p, const Real3& d) const
{
    if(this->_is_inside(p))
    {
        Real tmin = inf;
        for(std::size_t i=0; i<3; ++i)
        {
            if(std::abs(d[i]) < epsilon) continue;
            const Real tmp = (d[i] > 0) ? (this->upper_[i] - p[i]) / d[i] :
                                          (this->lower_[i] - p[i]) / d[i] ;
            tmin = std::min(tmin, tmp);
        }
        bool retval(tmin <= 1.0);
        return std::make_pair(retval, tmin);
    }
    else
    {
        Real tmin;
        Real3 q;
        const bool retval(collision::intersect_ray_AABB(p, d, lower_, upper_, tmin, q));
        return std::make_pair(retval, tmin);
    }
}

}// ecell4
