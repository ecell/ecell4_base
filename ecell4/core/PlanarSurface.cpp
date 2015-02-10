#include "PlanarSurface.hpp"
#include "collision.hpp"

namespace ecell4
{

PlanarSurface::PlanarSurface()
    : origin_(0,0,0), e0_(1,0,0), e1_(0,1,0), n_(0,0,1)
{
    ;
}

PlanarSurface::PlanarSurface(
        const Real3& origin, const Real3& e0, const Real3& e1)
    : origin_(origin), e0_(e0), e1_(e1)
{
    n_ = cross_product(e0_, e1_);
    n_ /= length(n_);
}

PlanarSurface::PlanarSurface(const PlanarSurface& rhs)
    : origin_(rhs.origin_), e0_(rhs.e0_), e1_(rhs.e1_), n_(rhs.n_)
{
    ;
}

Real PlanarSurface::is_inside(const Real3& coord) const
{
    return dot_product(origin_ - coord, n_);
}

Real3 PlanarSurface::draw_position(
    boost::shared_ptr<RandomNumberGenerator>& rng) const
{
    const Real a(rng->uniform(0,1)),
               b(rng->uniform(0,1));
    return origin_ + e0_ * a + e1_ * b;
}

bool PlanarSurface::test_AABB(const Real3& lower, const Real3& upper) const
{
    return collision::test_AABB_plane(AABB(lower, upper), *this);
}

}
