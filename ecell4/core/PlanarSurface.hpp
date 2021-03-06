#ifndef ECELL4_PLANAR_SURFACE_HPP
#define ECELL4_PLANAR_SURFACE_HPP

#include "Shape.hpp"

namespace ecell4
{

struct PlanarSurface
    : public Shape
{

    PlanarSurface();
    PlanarSurface(const Real3& origin, const Real3& e0, const Real3& e1);
    PlanarSurface(const PlanarSurface& rhs);
    Real is_inside(const Real3& coord) const;
    Real3 draw_position(
        std::shared_ptr<RandomNumberGenerator>& rng) const;
    bool test_AABB(const Real3& lower, const Real3& upper) const;
    void bounding_box(
        const Real3& edge_lengths, Real3& lower, Real3& u) const;

    const Real3& origin() const
    {
        return origin_;
    }

    const Real3& e0() const
    {
        return e0_;
    }

    const Real3& e1() const
    {
        return e1_;
    }

    const Real3& normal() const
    {
        return n_;
    }

    dimension_kind dimension() const
    {
        return TWO;
    }

protected:

    Real3 origin_, e0_, e1_, n_;
    Real d_;
};

}

#endif /* ECELL4_PLANAR_SURFACE_HPP */
