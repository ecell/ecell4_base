#ifndef __ECELL4_PLANAR_SURFACE_HPP
#define __ECELL4_PLANAR_SURFACE_HPP

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
        boost::shared_ptr<RandomNumberGenerator>& rng) const;

    dimension_kind dimension() const
    {
        return TWO;
    }

protected:

    Real3 origin_, e0_, e1_, n_;

};

}

#endif /* __ECELL4_PLANAR_SURFACE_HPP */
