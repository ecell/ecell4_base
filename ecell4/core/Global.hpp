#ifndef __ECELL4__GLOBAL_HPP
#define __ECELL4__GLOBAL_HPP

#include <ostream>
#include <iomanip>
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

inline bool operator<(const Global& lhs, const Global& rhs)
{
    return (lhs.col < rhs.col ? true :
        (lhs.row < rhs.row ? true : (lhs.layer < lhs.layer ? true : false)));
}

inline bool operator>(const Global& lhs, const Global& rhs)
{
    return (lhs.col > rhs.col ? true :
        (lhs.row > rhs.row ? true : (lhs.layer > lhs.layer ? true : false)));
}

inline bool operator==(const Global& lhs, const Global& rhs)
{
    return (lhs.col == rhs.col && lhs.row == rhs.row && lhs.layer == rhs.layer);
}

inline bool operator!=(const Global& lhs, const Global& rhs)
{
    return (lhs.col != rhs.col || lhs.row != rhs.row || lhs.layer != rhs.layer);
}

template<typename Tstrm_, typename Ttraits_>
inline std::basic_ostream<Tstrm_, Ttraits_>& operator<<(
    std::basic_ostream<Tstrm_, Ttraits_>& strm, const Global& g)
{
    strm << "{" << g.col <<  ", " << g.row <<  ", " << g.layer << "}";
    return strm;
}

} // ecell4

#endif
