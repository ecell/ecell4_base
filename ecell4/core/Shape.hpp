#ifndef __ECELL4_SHAPE_HPP
#define __ECELL4_SHAPE_HPP

#include "Position3.hpp"

namespace ecell4
{

struct Shape
{
    virtual bool func(const Position3& coord) = 0; // TODO: temporary name
};

} // ecell4

#endif /* __ECELL4_SHAPE_HPP */
