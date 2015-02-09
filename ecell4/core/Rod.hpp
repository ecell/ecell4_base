#ifndef __ECELL4_ROD_HPP
#define __ECELL4_ROD_HPP

#include "Shape.hpp"

namespace ecell4
{

struct RodSurface;

// A Rod is parallel with the x-axis.
// The center of a Rod is the origin.
struct Rod
    : public Shape
{
public:

    Rod();
    Rod(const Real& length, const Real& radius);
    Rod(const Rod& rhs);
    const Real& lengthX() const;
    const Real& radius() const;
    Real is_inside(const Real3& pos) const;
    Real distance(const Real3& pos) const;
    Real3 draw_position(boost::shared_ptr<RandomNumberGenerator>& rng) const;
    RodSurface surface() const;
    bool test_AABB(const Real3& l, const Real3& u) const;

protected:

    Real length_; // LengthX
    Real radius_; // LengthY/2
};

struct RodSurface
    : public Shape
{
public:

    RodSurface();
    RodSurface(const Real& length, const Real& radius);
    RodSurface(const RodSurface& rhs);
    const Real& lengthX() const;
    const Real& radius() const;
    Real is_inside(const Real3& pos) const;
    Real distance(const Real3& pos) const;
    Real3 draw_position(boost::shared_ptr<RandomNumberGenerator>& rng) const;
    Rod inside() const;
    bool test_AABB(const Real3& l, const Real3& u) const;

    dimension_kind dimension() const
    {
        return TWO;
    }

protected:

    Real length_; // LengthX
    Real radius_; // LengthY/2
};

} // ecell4

#endif /* __ECELL4_ROD_HPP */
