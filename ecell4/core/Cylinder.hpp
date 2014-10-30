#ifndef __ECELL4_CYLINDER_HPP
#define __ECELL4_CYLINDER_HPP

#include <ostream>
#include "Shape.hpp"

namespace ecell4
{

struct CylindricalSurface;

struct Cylinder
    : public Shape
{
public:

    /** for epdp
     */
    typedef Position3 position_type;
    typedef position_type::value_type length_type;
    typedef position_type::value_type value_type;

public:

    Cylinder();
    Cylinder(const Position3& center, const Real radius,
        const Position3& axis, const Real half_height);
    Cylinder(const Cylinder& rhs);

    const Real& radius() const;
    const Position3& center() const;
    const Real& half_height() const;
    const Position3& axis() const;

    Real is_inside(const Position3& coord) const;
    Real distance(const Position3& pos) const;
    CylindricalSurface surface() const;
    Position3 draw_position(
        boost::shared_ptr<RandomNumberGenerator>& rng) const;

    inline const Position3& position() const
    {
        return center();
    }

    Position3& position()
    {
        return center_;
    }

    inline const Real& size() const
    {
        return radius_;
    }

    Real& size()
    {
        return radius_;
    }

protected:

    inline std::pair<Real, Real> to_internal(const Position3& pos) const
    {
        const Position3 v(pos - center_);
        const Real z(dot_product(v, axis_));
        const Real r(length(v - multiply(axis_, z)));
        return std::make_pair(r, z);
    }

protected:

    Position3 center_;
    Real radius_;
    Position3 axis_;
    Real half_height_;
};

struct CylindricalSurface
    : public Shape
{
    CylindricalSurface();
    CylindricalSurface(const Position3& center, const Real radius,
        const Position3& axis, const Real half_height);
    CylindricalSurface(const CylindricalSurface& rhs);

    const Real& radius() const;
    const Position3& center() const;
    const Real& half_height() const;
    const Position3& axis() const;

    Real is_inside(const Position3& coord) const;
    Real distance(const Position3& pos) const;
    Cylinder inside() const;
    Position3 draw_position(
        boost::shared_ptr<RandomNumberGenerator>& rng) const;

    dimension_kind dimension() const
    {
        return TWO;
    }

protected:

    Position3 center_;
    Real radius_;
    Position3 axis_;
    Real half_height_;
};

} // ecell4

#endif /* __ECELL4_CYLINDER_HPP */

