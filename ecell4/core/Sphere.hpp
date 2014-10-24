#ifndef __ECELL4_SPHERE_HPP
#define __ECELL4_SPHERE_HPP

#include <ostream>
#include "Shape.hpp"

namespace ecell4
{

struct SphericalSurface;

struct Sphere
    : public Shape
{
public:

    /** for epdp
     */
    typedef Position3 position_type;
    typedef position_type::value_type length_type;
    typedef position_type::value_type value_type;

public:

    Sphere();
    Sphere(const Position3& center, const Real radius);
    Sphere(const Sphere& rhs);
    const Real& radius() const;
    const Position3& center() const;
    Real is_inside(const Position3& coord) const;
    Real distance(const Position3& pos) const;
    SphericalSurface surface() const;
    Position3 draw_position(
        boost::shared_ptr<RandomNumberGenerator>& rng) const;

    inline const Position3& position() const
    {
        return center();
    }

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
    const Real& radius() const;
    const Position3& center() const;
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
