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
    Real is_inside(const Position3& coord) const;
    Real distance(const Position3& pos) const;
    Position3 draw_position_inside(
        boost::shared_ptr<RandomNumberGenerator>& rng) const;

protected:

    Position3 center_;
    Real radius_;
};

} // ecell4

#endif /* __ECELL4_SPHERE_HPP */
