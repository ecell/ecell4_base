#ifndef __ECELL4_SPHERE_HPP
#define __ECELL4_SPHERE_HPP

#include "Shape.hpp"

namespace ecell4
{

struct Sphere : public Shape
{
    Sphere(const Position3& center, const Real radius);
    bool func(const Position3& coord); // TODO: temporary name

protected:
    Position3 center_;
    Real radius_;
};

} // ecell4

#endif /* __ECELL4_SPHERE_HPP */
