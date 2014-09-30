#include "Sphere.hpp"

namespace ecell4
{

Sphere::Sphere()
    : center_(), radius_()
{
    ;
}

Sphere::Sphere(const Position3& center, const Real radius)
    : center_(center), radius_(radius)
{
    ;
}

Sphere::Sphere(const Sphere& rhs)
    : center_(rhs.center()), radius_(rhs.radius())
{
    ;
}

Real Sphere::radius() const
{
    return radius_;
}

Position3 Sphere::center() const
{
    return center_;
}

Real Sphere::distance(const Position3& coord) const
{
    return length(coord - center_) - radius_;
}

Real Sphere::is_inside(const Position3& coord) const
{
    return distance(coord);
}

SphericalSurface Sphere::surface() const
{
    return SphericalSurface(center_, radius_);
}

Position3 Sphere::draw_position(
    boost::shared_ptr<RandomNumberGenerator>& rng) const
{
    if (radius_ <= 0.0)
    {
        return center_;
    }

    while (true)
    {
        const Real x(rng->uniform(-radius_, +radius_));
        const Real y(rng->uniform(-radius_, +radius_));
        const Real z(rng->uniform(-radius_, +radius_));
        const Position3 dir(x, y, z);
        const Position3 pos(dir + center_);
        if (is_inside(pos) <= 0.0)
        {
            return pos;
        }
    }

    ; // never reach here
}

SphericalSurface::SphericalSurface()
    : center_(), radius_()
{
    ;
}

SphericalSurface::SphericalSurface(const Position3& center, const Real radius)
    : center_(center), radius_(radius)
{
    ;
}

SphericalSurface::SphericalSurface(const SphericalSurface& rhs)
    : center_(rhs.center()), radius_(rhs.radius())
{
    ;
}

Real SphericalSurface::radius() const
{
    return radius_;
}

Position3 SphericalSurface::center() const
{
    return center_;
}

Real SphericalSurface::distance(const Position3& coord) const
{
    return length(coord - center_) - radius_;
}

Real SphericalSurface::is_inside(const Position3& coord) const
{
    return distance(coord);
}

Sphere SphericalSurface::inside() const
{
    return Sphere(center_, radius_);
}

Position3 SphericalSurface::draw_position(
    boost::shared_ptr<RandomNumberGenerator>& rng) const
{
    if (radius_ <= 0.0)
    {
        return center_;
    }

    return rng->direction3d(radius_) + center_;
}

} // ecell4
