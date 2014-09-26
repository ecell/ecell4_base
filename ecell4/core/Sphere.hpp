#ifndef __ECELL4_SPHERE_HPP
#define __ECELL4_SPHERE_HPP

#include "Shape.hpp"

namespace ecell4
{

struct SphericalSurface;

struct Sphere
    : public Shape
{
    Sphere();
    Sphere(const Position3& center, const Real radius);
    Sphere(const Sphere& rhs);
    Real radius() const;
    Position3 center() const;
    Real is_inside(const Position3& coord) const;
    Real distance(const Position3& pos) const;
    SphericalSurface surface() const;
    Position3 draw_position(
        boost::shared_ptr<RandomNumberGenerator>& rng) const;

protected:

    Position3 center_;
    Real radius_;
};

struct SphericalSurface
    : public Shape
{
    SphericalSurface();
    SphericalSurface(const Position3& center, const Real radius);
    SphericalSurface(const SphericalSurface& rhs);
    Real radius() const;
    Position3 center() const;
    Real is_inside(const Position3& coord) const;
    Real distance(const Position3& pos) const;
    Sphere inside() const;
    Position3 draw_position(
        boost::shared_ptr<RandomNumberGenerator>& rng) const;

    dimension_kind dimension() const
    {
        return TWO;
    }

protected:

    Position3 center_;
    Real radius_;
};

} // ecell4

#endif /* __ECELL4_SPHERE_HPP */
