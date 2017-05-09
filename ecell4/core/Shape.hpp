#ifndef ECELL4_SHAPE_HPP
#define ECELL4_SHAPE_HPP

#include "Real3.hpp"
#include "RandomNumberGenerator.hpp"


namespace ecell4
{

struct Shape
{
    typedef enum
    {
        ONE = 1,
        TWO = 2,
        THREE = 3,
        UNDEF = 4,
    } dimension_kind;

    virtual ~Shape()
    {
        ; // do nothing
    }

    virtual dimension_kind dimension() const = 0;
    // virtual dimension_kind dimension() const
    // {
    //     return THREE;
    // }

    virtual Real is_inside(const Real3& coord) const = 0;
    virtual Real3 draw_position(
        boost::shared_ptr<RandomNumberGenerator>& rng) const = 0;
    virtual bool test_AABB(const Real3& l, const Real3& u) const = 0;

    virtual void bounding_box(
        const Real3& edge_lengths, Real3& lower, Real3& upper) const
    {
        lower = Real3(0.0, 0.0, 0.0);
        upper = edge_lengths;
    }
};

} // ecell4

#endif /* ECELL4_SHAPE_HPP */
