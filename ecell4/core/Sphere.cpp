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

Real Sphere::distance(const Position3& coord) const
{
    return length(coord - center_) - radius_;
}

Real Sphere::is_inside(const Position3& coord) const
{
    return distance(coord);
}

Position3 Sphere::draw_position_inside(
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

} // ecell4
