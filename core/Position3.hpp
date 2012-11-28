#ifndef __POSITION3_HPP
#define __POSITION3_HPP

#include <ostream>
#include <iomanip>
#include <functional>
#include <algorithm>
#include <cmath>
#include <gsl/gsl_pow_int.h>
#include <boost/array.hpp>

#include "types.hpp"
#include "functions.hpp"


namespace ecell4
{

struct Position3
    : public boost::array<Real, 3>
{
    typedef boost::array<Real, 3> base_type;
    typedef typename base_type::value_type value_type;
    typedef typename base_type::size_type size_type;

    Position3& operator+=(Position3 const& rhs);
    Position3& operator-=(Position3 const& rhs);
    Position3& operator*=(Position3::value_type const& rhs);
    Position3& operator/=(Position3::value_type const& rhs);

    Position3()
    {
        (*this)[0] = 0;
        (*this)[1] = 0;
        (*this)[2] = 0;
    }

    Position3(value_type p0, value_type p1, value_type p2)
    {
        (*this)[0] = p0;
        (*this)[1] = p1;
        (*this)[2] = p2;
    }

    // Position3(const Real (&a)[3])
    //     : base_type(*reinterpret_cast<const base_type*>(&a))
    // {
    //     ;
    // }

    // Position3(const Real a[3])
    //     : base_type(*reinterpret_cast<const base_type*>(a))
    // {
    //     ;
    // }

    // Position3(const base_type& a)
    //     : base_type(a)
    // {
    //     ;
    // }
};

inline Position3 add(Position3 const& p1, Position3 const& p2)
{
    Position3 retval;
    retval[0] = p1[0] + p2[0];
    retval[1] = p1[1] + p2[1];
    retval[2] = p1[2] + p2[2];
    return retval;
}

inline Position3 subtract(Position3 const& p1, Position3 const& p2)
{
    Position3 retval;
    retval[0] = p1[0] - p2[0];
    retval[1] = p1[1] - p2[1];
    retval[2] = p1[2] - p2[2];
    return retval;
}

inline Position3 divide(Position3 const& p1, Position3::value_type const& p2)
{
    Position3 retval;
    retval[0] = p1[0] / p2;
    retval[1] = p1[1] / p2;
    retval[2] = p1[2] / p2;
    return retval;
}

inline Position3 multiply(Position3 const& p1, Position3::value_type const& p2)
{
    Position3 retval;
    retval[0] = p1[0] * p2;
    retval[1] = p1[1] * p2;
    retval[2] = p1[2] * p2;
    return retval;
}

inline Position3 modulo(Position3 const& p1, Position3::value_type const& p2)
{
    Position3 retval;
    retval[0] = modulo(p1[0], p2);
    retval[1] = modulo(p1[1], p2);
    retval[2] = modulo(p1[2], p2);
    return retval;
}

inline Position3 modulo(Position3 const& p1, Position3 const& p2)
{
    Position3 retval;
    retval[0] = modulo(p1[0], p2[0]);
    retval[1] = modulo(p1[1], p2[1]);
    retval[2] = modulo(p1[2], p2[2]);
    return retval;
}

inline Position3 abs(Position3 const& v)
{
    Position3 retval;
    retval[0] = abs(v[0]);
    retval[1] = abs(v[1]);
    retval[2] = abs(v[2]);
    return retval;
}

inline Position3::value_type dot_product(
    Position3 const& p1, Position3 const& p2)
{
    return p1[0] * p2[0] + p1[1] * p2[1] + p1[2] * p2[2];
}

inline Position3 cross_product(Position3 const& p1, Position3 const& p2)
{
    Position3 retval;
    retval[0] = p1[1] * p2[2] - p1[2] * p2[1];
    retval[1] = p1[2] * p2[0] - p1[0] * p2[2];
    retval[2] = p1[0] * p2[1] - p1[1] * p2[0];
    return retval;
}

inline Position3::size_type length_sq(Position3 const& r)
{
    return gsl_pow_2(r[0]) + gsl_pow_2(r[1]) + gsl_pow_2(r[2]);
}

inline Position3::size_type length(Position3 const& r)
{
    return std::sqrt(length_sq(r));
}

inline Position3 operator+(Position3 const& lhs, Position3 const& rhs)
{
    return add(lhs, rhs);
}

inline Position3 operator-(Position3 const& lhs, Position3 const& rhs)
{
    return subtract(lhs, rhs);
}

inline Position3 operator/(
    Position3 const& lhs, Position3::value_type const& rhs)
{
    return divide(lhs, rhs);
}

inline Position3 operator*(
    Position3 const& lhs, Position3::value_type const& rhs)
{
    return multiply(lhs, rhs);
}

template<typename Tstrm_, typename Ttraits_>
inline std::basic_ostream<Tstrm_, Ttraits_>& operator<<(
    std::basic_ostream<Tstrm_, Ttraits_>& strm, const Position3& v)
{
    strm << std::setprecision(12)
         << "(" << v[0] <<  ", " << v[1] <<  ", " << v[2] << ")";
    return strm;
}

} // ecell4

#endif /* __POSITION3_HPP */
