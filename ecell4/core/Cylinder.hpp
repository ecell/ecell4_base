#ifndef ECELL4_CYLINDER_HPP
#define ECELL4_CYLINDER_HPP

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
    typedef Real3 position_type;
    typedef position_type::value_type length_type;
    typedef position_type::value_type value_type;

public:

    Cylinder();
    Cylinder(const Real3& center, const Real radius,
        const Real3& axis, const Real half_height);
    Cylinder(const Cylinder& rhs);

    const Real& radius() const;
    const Real3& center() const;
    const Real& half_height() const;
    const Real3& axis() const;

    Real is_inside(const Real3& coord) const;
    Real distance(const Real3& pos) const;
    bool test_AABB(const Real3& l, const Real3& u) const;
    CylindricalSurface surface() const;
    Real3 draw_position(
        boost::shared_ptr<RandomNumberGenerator>& rng) const;

    inline const Real3& position() const
    {
        return center();
    }

    Real3& position()
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

    inline std::pair<Real, Real> to_internal(const Real3& pos) const
    {
        const Real3 v(pos - center_);
        const Real z(dot_product(v, axis_));
        const Real r(length(v - multiply(axis_, z)));
        return std::make_pair(r, z);
    }

    dimension_kind dimension() const
    {
        return THREE;
    }

protected:

    Real3 center_;
    Real radius_;
    Real3 axis_;
    Real half_height_;
};

struct CylindricalSurface
    : public Shape
{
    CylindricalSurface();
    CylindricalSurface(const Real3& center, const Real radius,
        const Real3& axis, const Real half_height);
    CylindricalSurface(const CylindricalSurface& rhs);

    const Real& radius() const;
    const Real3& center() const;
    const Real& half_height() const;
    const Real3& axis() const;

    Real is_inside(const Real3& coord) const;
    Real distance(const Real3& pos) const;
    Cylinder inside() const;
    Real3 draw_position(
        boost::shared_ptr<RandomNumberGenerator>& rng) const;
    bool test_AABB(const Real3& l, const Real3& u) const;

    dimension_kind dimension() const
    {
        return TWO;
    }

protected:

    Real3 center_;
    Real radius_;
    Real3 axis_;
    Real half_height_;
};

} // ecell4

#endif /* ECELL4_CYLINDER_HPP */

