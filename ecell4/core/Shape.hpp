#ifndef __ECELL4_SHAPE_HPP
#define __ECELL4_SHAPE_HPP

#include "Real3.hpp"
#include "RandomNumberGenerator.hpp"


namespace ecell4
{

struct Shape
{
    typedef enum {
        UNDEF, ONE, TWO, THREE,
    } dimension_kind;

    virtual ~Shape()
    {
        ; // do nothing
    }

    virtual dimension_kind dimension() const
    {
        return THREE;
    }

    virtual Real is_inside(const Real3& coord) const = 0;
    virtual Real3 draw_position(
        boost::shared_ptr<RandomNumberGenerator>& rng) const = 0;
};

} // ecell4

#endif /* __ECELL4_SHAPE_HPP */
