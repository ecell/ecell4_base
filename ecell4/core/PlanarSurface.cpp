#include "PlanarSurface.hpp"

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
    const Real3 c(multiply(lower + upper, 0.5));
    const Real3 e(upper - c);

    const Real3 m(cross_product(e0_, e1_));
    const Real3 n(divide(m, length(m)));
    const Real d(dot_product(origin_, n));

    const Real r(e[0] * abs(n[0]) + e[1] * abs(n[1]) + e[2] * abs(n[2]));
    const Real s(dot_product(n, c) - d);
    return (abs(s) <= r);
}

}
