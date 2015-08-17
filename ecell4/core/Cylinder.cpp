#include "Cylinder.hpp"
#include "exceptions.hpp"
#include "collision.hpp"


namespace ecell4
{

Cylinder::Cylinder()
    : center_(), radius_(), axis_(), half_height_()
{
    ;
}

Cylinder::Cylinder(const Real3& center, const Real radius,
    const Real3& axis, const Real half_height)
    : center_(center), radius_(radius), axis_(axis), half_height_(half_height)
{
    ;
}

Cylinder::Cylinder(const Cylinder& rhs)
    : center_(rhs.center()), radius_(rhs.radius()),
    axis_(rhs.axis()), half_height_(rhs.half_height_)
{
    ;
}

const Real& Cylinder::radius() const
{
    return radius_;
}

const Real3& Cylinder::center() const
{
    return center_;
}

const Real& Cylinder::half_height() const
{
    return half_height_;
}

const Real3& Cylinder::axis() const
{
    return axis_;
}

Real Cylinder::distance(const Real3& coord) const
{
    return collision::distance_point_cylinder(coord, *this);
}

Real Cylinder::is_inside(const Real3& coord) const
{
    return distance(coord);
}

CylindricalSurface Cylinder::surface() const
{
    return CylindricalSurface(center_, radius_, axis_, half_height_);
}

Real3 Cylinder::draw_position(
    boost::shared_ptr<RandomNumberGenerator>& rng) const
{
    throw NotImplemented("not implemented yet.");
}

bool Cylinder::test_AABB(const Real3& l, const Real3& u) const
{
    throw NotImplemented("not implemented yet.");
}

CylindricalSurface::CylindricalSurface()
    : center_(), radius_(), axis_(), half_height_()
{
    ;
}

CylindricalSurface::CylindricalSurface(const Real3& center, const Real radius,
    const Real3& axis, const Real half_height)
    : center_(center), radius_(radius), axis_(axis), half_height_(half_height_)
{
    ;
}

CylindricalSurface::CylindricalSurface(const CylindricalSurface& rhs)
    : center_(rhs.center()), radius_(rhs.radius()),
    axis_(rhs.axis()), half_height_(rhs.half_height_)
{
    ;
}

const Real& CylindricalSurface::radius() const
{
    return radius_;
}

const Real3& CylindricalSurface::center() const
{
    return center_;
}

const Real& CylindricalSurface::half_height() const
{
    return half_height_;
}

const Real3& CylindricalSurface::axis() const
{
    return axis_;
}

Real CylindricalSurface::distance(const Real3& coord) const
{
    return inside().distance(coord); //XXX: This is too slow.
}

Real CylindricalSurface::is_inside(const Real3& coord) const
{
    return distance(coord);
}

Cylinder CylindricalSurface::inside() const
{
    return Cylinder(center_, radius_, axis_, half_height_);
}

Real3 CylindricalSurface::draw_position(
    boost::shared_ptr<RandomNumberGenerator>& rng) const
{
    throw NotImplemented("not implemented yet.");
}

bool CylindricalSurface::test_AABB(const Real3& l, const Real3& u) const
{
    throw NotImplemented("not implemented yet.");
}

} // ecell4

