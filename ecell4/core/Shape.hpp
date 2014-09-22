#ifndef __ECELL4_SHAPE_HPP
#define __ECELL4_SHAPE_HPP

#include "Position3.hpp"

namespace ecell4
{

struct Shape
{
    virtual ~Shape()
    {
        ; // do nothing
    }

    virtual bool is_inside(const Position3& coord) const = 0;
};

} // ecell4

#endif /* __ECELL4_SHAPE_HPP */
