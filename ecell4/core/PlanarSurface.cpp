#include "PlanarSurface.hpp"
#include "collision.hpp"

namespace ecell4
{

PlanarSurface::PlanarSurface()
    : origin_(0, 0, 0), e0_(1, 0, 0), e1_(0, 1, 0), n_(0, 0, 1), d_(0.0)
{
    ;
}

PlanarSurface::PlanarSurface(
        const Real3& origin, const Real3& e0, const Real3& e1)
    : origin_(origin), e0_(e0), e1_(e1)
{
    n_ = cross_product(e0_, e1_);
    n_ /= length(n_);
    d_ = dot_product(origin_, n_);
}

PlanarSurface::PlanarSurface(const PlanarSurface& rhs)
    : origin_(rhs.origin_), e0_(rhs.e0_), e1_(rhs.e1_), n_(rhs.n_), d_(rhs.d_)
{
    ;
}

Real PlanarSurface::is_inside(const Real3& coord) const
{
    return d_ - dot_product(coord, n_);
}

Real3 PlanarSurface::draw_position(
    boost::shared_ptr<RandomNumberGenerator>& rng) const
{
    const Real a(rng->uniform(0, 1)),
               b(rng->uniform(0, 1));
    return origin_ + e0_ * a + e1_ * b;
}

bool PlanarSurface::test_AABB(const Real3& lower, const Real3& upper) const
{
    return collision::test_AABB_plane(AABB(lower, upper), *this);
}

void PlanarSurface::bounding_box(
    const Real3& edge_lengths, Real3& lower, Real3& upper) const
{
    if (n_[0] > epsilon)
    {
        lower[0] = std::max(
            (d_ - n_[1] * edge_lengths[1] - n_[2] * edge_lengths[2]) / n_[0], 0.0);
        upper[0] = std::min(d_ / n_[0], edge_lengths[0]);
    }
    else
    {
        lower[0] = 0.0;
        upper[0] = edge_lengths[0];
    }
    if (n_[1] > epsilon)
    {
        lower[1] = std::max(
            (d_ - n_[0] * edge_lengths[0] - n_[2] * edge_lengths[2]) / n_[1], 0.0);
        upper[1] = std::min(d_ / n_[1], edge_lengths[1]);
    }
    else
    {
        lower[1] = 0.0;
        upper[1] = edge_lengths[1];
    }
    if (n_[2] > epsilon)
    {
        lower[2] = std::max(
            (d_ - n_[1] * edge_lengths[1] - n_[0] * edge_lengths[0]) / n_[2], 0.0);
        upper[2] = std::min(d_ / n_[2], edge_lengths[2]);
    }
    else
    {
        lower[2] = 0.0;
        upper[2] = edge_lengths[2];
    }
}

}
