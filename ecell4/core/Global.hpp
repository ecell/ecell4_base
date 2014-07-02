#ifndef __ECELL4__GLOBAL_HPP
#define __ECELL4__GLOBAL_HPP

#include <vector>
#include "types.hpp"

namespace ecell4
{

struct Global
{
    Integer col;
    Integer row;
    Integer layer;

    Global()
    {
        this->col = 0;
        this->row = 0;
        this->layer = 0;
    }

    Global(Integer col, Integer row, Integer layer)
    {
        this->col = col;
        this->row = row;
        this->layer = layer;
    }

    Global(const Global& global)
    {
        this->col = global.col;
        this->row = global.row;
        this->layer = global.layer;
    }

    Global east() const;
    Global west() const;
    Global south() const;
    Global north() const;
    Global dorsal() const;
    Global ventral() const;

};

inline Global add(const Global& g1, const Global& g2)
{
    Global retval;
    retval.col = g1.col + g2.col;
    retval.row = g1.row + g2.row;
    retval.layer = g1.layer + g2.layer;
    return retval;
}

inline Global subtract(const Global& g1, const Global& g2)
{
    Global retval;
    retval.col = g1.col - g2.col;
    retval.row = g1.row - g2.row;
    retval.layer = g1.layer - g2.layer;
    return retval;
}

inline Global operator+(const Global& lhs, const Global& rhs)
{
    return add(lhs, rhs);
}

inline Global operator-(const Global& lhs, const Global& rhs)
{
    return subtract(lhs, rhs);
}

} // ecell4

#endif
