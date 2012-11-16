#ifndef __VECTOR3_HPP
#define __VECTOR3_HPP

#include <ostream>
#include <iomanip>
#include <functional>
#include <algorithm>
#include <boost/array.hpp>

#include "types.hpp"
// #include "linear_algebra.hpp"


namespace ecell4
{

struct Position3
    : public boost::array<Real, 3>
{
    typedef boost::array<Real, 3> base_type;
    typedef typename base_type::value_type value_type;
    typedef typename base_type::size_type size_type;

    // Position3& operator+=(Position3 const& rhs)
    // {
    //     *this = add(*this, rhs);
    //     return *this;
    // }

    // Position3& operator-=(Position3 const& rhs)
    // {
    //     *this = subtract(*this, rhs);
    //     return *this;
    // }

    // template<typename TT_>
    // Position3& operator*=(TT_ const& rhs)
    // {
    //     *this = multiply(*this, rhs);
    //     return *this;
    // }

    // template<typename TT_>
    // Position3& operator/=(TT_ const& rhs)
    // {
    //     *this = divide(*this, rhs);
    //     return *this;
    // }

    Position3()
    {
        (*this)[0] = 0;
        (*this)[1] = 0;
        (*this)[2] = 0;
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

    Position3(value_type p0, value_type p1, value_type p2)
    {
        (*this)[0] = p0;
        (*this)[1] = p1;
        (*this)[2] = p2;
    }
};

// inline Position3 operator+(Position3 const& lhs, Position3 const& rhs)
// {
//     return add(lhs, rhs);
// }

// inline Position3 operator-(Position3 const& lhs, Position3 const& rhs)
// {
//     return subtract(lhs, rhs);
// }

// inline Position3 operator/(
//     Position3 const& lhs, Position3::value_type const& rhs)
// {
//     return divide(lhs, rhs);
// }

// inline Position3 operator*(
//     Position3 const& lhs, Position3::value_type const& rhs)
// {
//     return multiply(lhs, rhs);
// }

template<typename Tstrm_, typename Ttraits_>
inline std::basic_ostream<Tstrm_, Ttraits_>& operator<<(
    std::basic_ostream<Tstrm_, Ttraits_>& strm, const Position3& v)
{
    strm << std::setprecision(12)
         << "(" << v[0] <<  ", " << v[1] <<  ", " << v[2] << ")";
    return strm;
}

} // ecell4

#endif /* __VECTOR3_HPP */
