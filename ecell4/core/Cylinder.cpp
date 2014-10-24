#include "Cylinder.hpp"
#include "exceptions.hpp"


namespace ecell4
{

Cylinder::Cylinder()
    : center_(), radius_(), axis_(), half_height_()
{
    ;
}

Cylinder::Cylinder(const Position3& center, const Real radius,
    const Position3& axis, const Real half_height)
    : center_(center), radius_(radius), axis_(axis), half_height_(half_height_)
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

const Position3& Cylinder::center() const
{
    return center_;
}

const Real& Cylinder::half_height() const
{
    return half_height_;
}

const Position3& Cylinder::axis() const
{
    return axis_;
}

Real Cylinder::distance(const Position3& coord) const
{
    /* First compute the (z,r) components of pos in a coordinate system 
     * defined by the vectors unitR and unit_z, where unitR is
     * choosen such that unitR and unit_z define a plane in which
     * pos lies. */
    const std::pair<Real, Real> r_z(to_internal(coord));

    /* Then compute distance to cylinder. */
    const Real dz(std::fabs(r_z.second) - half_height_);
    const Real dr(r_z.first - radius_);

    Real L;
    if (dz > 0)
    {
        // pos is (either) to the right or to the left of the cylinder.
        if (r_z.first > radius_)
        {
            // Compute distance to edge.
            L = std::sqrt(dz * dz + dr * dr);
        }
        else
        {
            L = dz;
        }
    }
    else
    {
        if (dr > radius_)
        {
            // pos is somewhere 'parallel' to the cylinder.
            L = dr;
        }
        else
        {
            // Inside cylinder.
            L = std::max(dr, dz);
        }
    }
    return L;
}

Real Cylinder::is_inside(const Position3& coord) const
{
    return distance(coord);
}

CylindricalSurface Cylinder::surface() const
{
    return CylindricalSurface(center_, radius_, axis_, half_height_);
}

Position3 Cylinder::draw_position(
    boost::shared_ptr<RandomNumberGenerator>& rng) const
{
    throw NotImplemented("not implemented yet.");
}

CylindricalSurface::CylindricalSurface()
    : center_(), radius_(), axis_(), half_height_()
{
    ;
}

CylindricalSurface::CylindricalSurface(const Position3& center, const Real radius,
    const Position3& axis, const Real half_height)
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

const Position3& CylindricalSurface::center() const
{
    return center_;
}

const Real& CylindricalSurface::half_height() const
{
    return half_height_;
}

const Position3& CylindricalSurface::axis() const
{
    return axis_;
}

Real CylindricalSurface::distance(const Position3& coord) const
{
    return inside().distance(coord); //XXX: This is too slow.
}

Real CylindricalSurface::is_inside(const Position3& coord) const
{
    return distance(coord);
}

Cylinder CylindricalSurface::inside() const
{
    return Cylinder(center_, radius_, axis_, half_height_);
}

Position3 CylindricalSurface::draw_position(
    boost::shared_ptr<RandomNumberGenerator>& rng) const
{
    throw NotImplemented("not implemented yet.");
}

} // ecell4

