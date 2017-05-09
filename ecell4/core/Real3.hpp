#ifndef ECELL4_REAL3_HPP
#define ECELL4_REAL3_HPP

#include <ostream>
#include <iomanip>
#include <functional>
#include <algorithm>
#include <cmath>
#include <boost/array.hpp>

#include <ecell4/core/config.h>

#include "types.hpp"
#include "functions.hpp"

#include "hash.hpp"

namespace ecell4
{

struct Real3
    : public boost::array<Real, 3>
{
    typedef boost::array<Real, 3> base_type;
    typedef base_type::value_type value_type;
    typedef base_type::size_type size_type;

    Real3& operator+=(const Real3& rhs);
    Real3& operator-=(const Real3& rhs);
    Real3& operator*=(const Real3::value_type& rhs);
    Real3& operator/=(const Real3::value_type& rhs);

    Real3()
    {
        (*this)[0] = 0;
        (*this)[1] = 0;
        (*this)[2] = 0;
    }

    Real3(value_type p0, value_type p1, value_type p2)
    {
        (*this)[0] = p0;
        (*this)[1] = p1;
        (*this)[2] = p2;
    }

    Real3(const Real3 &rhs)
    {
        (*this)[0] = rhs[0];
        (*this)[1] = rhs[1];
        (*this)[2] = rhs[2];
    }

    // Real3(const Real (&a)[3])
    //     : base_type(*reinterpret_cast<const base_type*>(&a))
    // {
    //     ;
    // }

    // Real3(const Real a[3])
    //     : base_type(*reinterpret_cast<const base_type*>(a))
    // {
    //     ;
    // }

    // Real3(const base_type& a)
    //     : base_type(a)
    // {
    //     ;
    // }
};

inline Real3 add(const Real3& p1, const Real3& p2)
{
    Real3 retval;
    retval[0] = p1[0] + p2[0];
    retval[1] = p1[1] + p2[1];
    retval[2] = p1[2] + p2[2];
    return retval;
}

inline Real3 subtract(const Real3& p1, const Real3& p2)
{
    Real3 retval;
    retval[0] = p1[0] - p2[0];
    retval[1] = p1[1] - p2[1];
    retval[2] = p1[2] - p2[2];
    return retval;
}

inline Real3 divide(const Real3& p1, const Real3::value_type& p2)
{
    Real3 retval;
    retval[0] = p1[0] / p2;
    retval[1] = p1[1] / p2;
    retval[2] = p1[2] / p2;
    return retval;
}

inline Real3 multiply(const Real3& p1, const Real3::value_type& p2)
{
    Real3 retval;
    retval[0] = p1[0] * p2;
    retval[1] = p1[1] * p2;
    retval[2] = p1[2] * p2;
    return retval;
}

inline Real3 modulo(const Real3& p1, const Real3::value_type& p2)
{
    Real3 retval;
    retval[0] = modulo(p1[0], p2);
    retval[1] = modulo(p1[1], p2);
    retval[2] = modulo(p1[2], p2);
    return retval;
}

inline Real3 modulo(const Real3& p1, const Real3& p2)
{
    Real3 retval;
    retval[0] = modulo(p1[0], p2[0]);
    retval[1] = modulo(p1[1], p2[1]);
    retval[2] = modulo(p1[2], p2[2]);
    return retval;
}

inline Real3 abs(const Real3& v)
{
    Real3 retval;
    retval[0] = abs(v[0]);
    retval[1] = abs(v[1]);
    retval[2] = abs(v[2]);
    return retval;
}

inline Real3::value_type dot_product(
    const Real3& p1, const Real3& p2)
{
    return p1[0] * p2[0] + p1[1] * p2[1] + p1[2] * p2[2];
}

inline Real3 cross_product(const Real3& p1, const Real3& p2)
{
    Real3 retval;
    retval[0] = p1[1] * p2[2] - p1[2] * p2[1];
    retval[1] = p1[2] * p2[0] - p1[0] * p2[2];
    retval[2] = p1[0] * p2[1] - p1[1] * p2[0];
    return retval;
}

inline Real3::value_type length_sq(const Real3& r)
{
    return pow_2(r[0]) + pow_2(r[1]) + pow_2(r[2]);
}

inline Real3::value_type length(const Real3& r)
{
    return std::sqrt(length_sq(r));
}

inline Real3 operator+(const Real3& lhs, const Real3& rhs)
{
    return add(lhs, rhs);
}

inline Real3 operator-(const Real3& lhs, const Real3& rhs)
{
    return subtract(lhs, rhs);
}

inline Real3 operator/(
    const Real3& lhs, const Real3::value_type& rhs)
{
    return divide(lhs, rhs);
}

inline Real3 operator*(
    const Real3& lhs, const Real3::value_type& rhs)
{
    return multiply(lhs, rhs);
}

template<typename Tstrm_, typename Ttraits_>
inline std::basic_ostream<Tstrm_, Ttraits_>& operator<<(
    std::basic_ostream<Tstrm_, Ttraits_>& strm, const Real3& v)
{
    strm << std::setprecision(12)
         << "(" << v[0] <<  ", " << v[1] <<  ", " << v[2] << ")";
    return strm;
}

inline Real3 ones()
{
    return Real3(1.0, 1.0, 1.0);
}

inline Real3 unitx()
{
    return Real3(1.0, 0.0, 0.0);
}

inline Real3 unity()
{
    return Real3(0.0, 1.0, 0.0);
}

inline Real3 unitz()
{
    return Real3(0.0, 0.0, 1.0);
}

} // ecell4

ECELL4_DEFINE_HASH_BEGIN()

template<>
struct hash<ecell4::Real3>
{
    typedef ecell4::Real3 argument_type;

    std::size_t operator()(argument_type const& val)
    {
        return hash<argument_type::value_type>()(val[0]) ^
            hash<argument_type::value_type>()(val[1]) ^
            hash<argument_type::value_type>()(val[2]);
    }
};

ECELL4_DEFINE_HASH_END()

#endif /* ECELL4_REAL3_HPP */
