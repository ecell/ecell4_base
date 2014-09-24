#ifndef __ECELL4_SHAPE_HPP
#define __ECELL4_SHAPE_HPP

#include "Position3.hpp"
#include "RandomNumberGenerator.hpp"


namespace ecell4
{

struct Shape
{
    virtual ~Shape()
    {
        ; // do nothing
    }

    virtual Real is_inside(const Position3& coord) const = 0;
    virtual Position3 draw_position_inside(
        boost::shared_ptr<RandomNumberGenerator>& rng) const = 0;
};

} // ecell4

#endif /* __ECELL4_SHAPE_HPP */
