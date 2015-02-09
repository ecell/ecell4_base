#include "Sphere.hpp"
#include "AABB.hpp"


namespace ecell4
{

Sphere::Sphere()
    : center_(), radius_()
{
    ;
}

Sphere::Sphere(const Real3& center, const Real radius)
    : center_(center), radius_(radius)
{
    ;
}

Sphere::Sphere(const Sphere& rhs)
    : center_(rhs.center()), radius_(rhs.radius())
{
    ;
}

const Real& Sphere::radius() const
{
    return radius_;
}

const Real3& Sphere::center() const
{
    return center_;
}

Real Sphere::distance(const Real3& coord) const
{
    return length(coord - center_) - radius_;
}

Real Sphere::is_inside(const Real3& coord) const
{
    return distance(coord);
}

SphericalSurface Sphere::surface() const
{
    return SphericalSurface(center_, radius_);
}

Real3 Sphere::draw_position(
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
        const Real3 dir(x, y, z);
        const Real3 pos(dir + center_);
        if (is_inside(pos) <= 0.0)
        {
            return pos;
        }
    }

    ; // never reach here
}

bool Sphere::test_AABB(const Real3& l, const Real3& u) const
{
    const AABB b(l, u);
    const Real Lsq(b.distance_sq(center_));
    return (Lsq <= radius_ * radius_);
}

SphericalSurface::SphericalSurface()
    : center_(), radius_()
{
    ;
}

SphericalSurface::SphericalSurface(const Real3& center, const Real radius)
    : center_(center), radius_(radius)
{
    ;
}

SphericalSurface::SphericalSurface(const SphericalSurface& rhs)
    : center_(rhs.center()), radius_(rhs.radius())
{
    ;
}

const Real& SphericalSurface::radius() const
{
    return radius_;
}

const Real3& SphericalSurface::center() const
{
    return center_;
}

Real SphericalSurface::distance(const Real3& coord) const
{
    return length(coord - center_) - radius_;
}

Real SphericalSurface::is_inside(const Real3& coord) const
{
    return distance(coord);
}

Sphere SphericalSurface::inside() const
{
    return Sphere(center_, radius_);
}

Real3 SphericalSurface::draw_position(
    boost::shared_ptr<RandomNumberGenerator>& rng) const
{
    if (radius_ <= 0.0)
    {
        return center_;
    }

    return rng->direction3d(radius_) + center_;
}

bool SphericalSurface::test_AABB(const Real3& l, const Real3& u) const
{
    const AABB b(l, u);
    const Real rsq(radius_ * radius_);
    if (b.distance_sq(center_) > rsq)
    {
        return false;
    }
    else if (b.farthest_distance_sq(center_) < rsq)
    {
        return false;
    }
    return true;
}

} // ecell4
