#include "Sphere.hpp"

namespace ecell4
{

Sphere::Sphere(const Position3& center, const Real radius)
    : center_(center), radius_(radius)
{
    ;
}

bool Sphere::func(const Position3& coord) // TODO: temporary name
{
    return length(coord-center_) <= radius_;
}

} // ecell4
