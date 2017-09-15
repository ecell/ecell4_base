#ifndef ECELL4_SPHERE_HPP
#define ECELL4_SPHERE_HPP

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
    typedef Real3 position_type;
    typedef position_type::value_type length_type;
    typedef position_type::value_type value_type;

public:

    Sphere();
    Sphere(const Real3& center, const Real radius);
    Sphere(const Sphere& rhs);
    const Real& radius() const;
    const Real3& center() const;
    Real is_inside(const Real3& coord) const;
    Real distance(const Real3& pos) const;
    SphericalSurface surface() const;
    Real3 draw_position(
        boost::shared_ptr<RandomNumberGenerator>& rng) const;
    bool test_AABB(const Real3& l, const Real3& u) const;

    inline const Real3& position() const
    {
        return center_;
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

    dimension_kind dimension() const
    {
        return THREE;
    }

protected:

    Real3 center_;
    Real radius_;
};

struct SphericalSurface
    : public Shape
{
    SphericalSurface();
    SphericalSurface(const Real3& center, const Real radius);
    SphericalSurface(const SphericalSurface& rhs);
    const Real& radius() const;
    const Real3& center() const;
    Real is_inside(const Real3& coord) const;
    Real distance(const Real3& pos) const;
    Sphere inside() const;
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
};

} // ecell4

#endif /* ECELL4_SPHERE_HPP */
