#ifndef ECELL4_AABB_SURFACE_HPP
#define ECELL4_AABB_SURFACE_HPP

#include "Shape.hpp"
#include "shape_operators.hpp"

namespace ecell4
{

struct AABBSurface
    : public Shape
{
    AABBSurface()
        : lower_(), upper_()
    {
        ;
    }

    AABBSurface(const Real3& lower, const Real3& upper)
        : lower_(lower), upper_(upper)
    {
        ;
    }

    AABBSurface(const AABBSurface& rhs)
        : lower_(rhs.lower()), upper_(rhs.upper())
    {
        ;
    }

    const Real3& lower() const
    {
        return lower_;
    }

    const Real3& upper() const
    {
        return upper_;
    }

    const Real3 center() const
    {
        return multiply(upper_ + lower_, 0.5);
    }

    const Real3 radius() const
    {
        return multiply(upper_ - lower_, 0.5);
    }

    Real distance_sq(const Real3 pos) const;
    Real distance(const Real3& pos) const;

    Real is_inside(const Real3& coord) const
    {
        if(this->_is_inside(coord))
            return -1. * this->distance(coord);
        else
            return this->distance(coord);
    }

    bool _is_inside(const Real3& coord) const
    {
        return (lower_[0] <= coord[0] && coord[0] <= upper_[0]) &&
               (lower_[1] <= coord[1] && coord[1] <= upper_[1]) &&
               (lower_[2] <= coord[2] && coord[2] <= upper_[2]);
    }

    Real3 draw_position(boost::shared_ptr<RandomNumberGenerator>& rng) const;
    bool test_AABB(const Real3& l, const Real3& u) const;
    bool test_segment(const Real3& p0, const Real3& p1) const;
    std::pair<bool, Real> intersect_ray(const Real3& p, const Real3& d) const;

    bool test_ray(const Real3& p, const Real3& d) const
    {
        return intersect_ray(p, d).first;
    }

    inline Real3 corner(const int& n) const
    {
        const Real3 p(
            ((n & 1) ? upper_[0] : lower_[0]),
            ((n & 2) ? upper_[1] : lower_[1]),
            ((n & 4) ? upper_[2] : lower_[2]));
        return p;
    }

    dimension_kind dimension() const
    {
        return TWO;
    }

    Surface surface() const
    {
        return Surface(boost::shared_ptr<Shape>(new AABBSurface(*this)));
    }

protected:

    Real3 lower_, upper_;
};

}// ecell4

#endif /* ECELL4_AABB_SURFACE_HPP */
