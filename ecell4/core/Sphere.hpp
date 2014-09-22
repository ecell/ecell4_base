#ifndef __ECELL4_SPHERE_HPP
#define __ECELL4_SPHERE_HPP

#include "Shape.hpp"

namespace ecell4
{

struct Sphere
    : public Shape
{
    Sphere();
    Sphere(const Position3& center, const Real radius);
    bool is_inside(const Position3& coord) const;

protected:

    Position3 center_;
    Real radius_;
};

} // ecell4

#endif /* __ECELL4_SPHERE_HPP */
