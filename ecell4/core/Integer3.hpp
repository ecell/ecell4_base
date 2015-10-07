#ifndef __ECELL4__GLOBAL_HPP
#define __ECELL4__GLOBAL_HPP

#include <ostream>
#include <iomanip>
#include <cmath>
#include <vector>

#include "types.hpp"
#include "exceptions.hpp"
#include "Real3.hpp"


namespace ecell4
{

struct Integer3
{
    typedef Integer value_type;
    typedef std::size_t size_type;

    value_type col;
    value_type row;
    value_type layer;

    Integer3()
    {
        this->col = 0;
        this->row = 0;
        this->layer = 0;
    }

    Integer3(value_type col, value_type row, value_type layer)
    {
        this->col = col;
        this->row = row;
        this->layer = layer;
    }

    Integer3(const Integer3& global)
    {
        this->col = global.col;
        this->row = global.row;
        this->layer = global.layer;
    }

    Integer3 east() const;
    Integer3 west() const;
    Integer3 south() const;
    Integer3 north() const;
    Integer3 dorsal() const;
    Integer3 ventral() const;

    value_type& operator[](size_type i)
    {
        switch (i)
        {
        case 0:
            return this->col;
        case 1:
            return this->row;
        case 2:
            return this->layer;
        }
        throw NotSupported("out of range");
    }

    const value_type& operator[](size_type i) const
    {
        switch (i)
        {
        case 0:
            return this->col;
        case 1:
            return this->row;
        case 2:
            return this->layer;
        }
        throw NotSupported("out of range");
    }
};

inline Integer3 add(const Integer3& g1, const Integer3& g2)
{
    Integer3 retval;
    retval.col = g1.col + g2.col;
    retval.row = g1.row + g2.row;
    retval.layer = g1.layer + g2.layer;
    return retval;
}

inline Integer3 subtract(const Integer3& g1, const Integer3& g2)
{
    Integer3 retval;
    retval.col = g1.col - g2.col;
    retval.row = g1.row - g2.row;
    retval.layer = g1.layer - g2.layer;
    return retval;
}

inline Integer3 abs(const Integer3& g1)
{
    Integer3 retval;
    retval.col = std::abs(g1.col);
    retval.row = std::abs(g1.row);
    retval.layer = std::abs(g1.layer);
    return retval;
}

inline Integer3 operator+(const Integer3& lhs, const Integer3& rhs)
{
    return add(lhs, rhs);
}

inline Integer3 operator-(const Integer3& lhs, const Integer3& rhs)
{
    return subtract(lhs, rhs);
}

inline bool operator<(const Integer3& lhs, const Integer3& rhs)
{
    return (lhs.col < rhs.col ? true :
        (lhs.row < rhs.row ? true : (lhs.layer < lhs.layer ? true : false)));
}

inline bool operator>(const Integer3& lhs, const Integer3& rhs)
{
    return (lhs.col > rhs.col ? true :
        (lhs.row > rhs.row ? true : (lhs.layer > lhs.layer ? true : false)));
}

inline bool operator==(const Integer3& lhs, const Integer3& rhs)
{
    return (lhs.col == rhs.col && lhs.row == rhs.row && lhs.layer == rhs.layer);
}

inline bool operator!=(const Integer3& lhs, const Integer3& rhs)
{
    return (lhs.col != rhs.col || lhs.row != rhs.row || lhs.layer != rhs.layer);
}

template<typename Tstrm_, typename Ttraits_>
inline std::basic_ostream<Tstrm_, Ttraits_>& operator<<(
    std::basic_ostream<Tstrm_, Ttraits_>& strm, const Integer3& g)
{
    strm << "{" << g.col <<  ", " << g.row <<  ", " << g.layer << "}";
    return strm;
}

} // ecell4

#endif
