#ifndef __LINEAR_ALGEBRA_HPP
#define __LINEAR_ALGEBRA_HPP

#include <algorithm>
#include <cmath>
#include <gsl/gsl_pow_int.h>

#include "Position3.hpp"


namespace ecell4
{

// inline Integer modulo(Integer const& p1, Integer const& p2)
// {
//     Integer r = p1 % p2;
//     if (r != 0 && (r > 0) == (p2 < 0))
//         r += p2;
//     return r;
// }

inline Real modulo(Real const& p1, Real const& p2)
{
    Real r = std::fmod(p1, p2);
    if (r != 0 && (r > 0) == (p2 < 0))
        r += p2;
    return r;
}

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
    retval[0] = multiply(p1[0], p2);
    retval[1] = multiply(p1[1], p2);
    retval[2] = multiply(p1[2], p2);
    return retval;
}

inline Position3 modulo(Position3 const& p1, Position3::value_type const& p2)
{
    T_ retval;
    retval[0] = modulo(p1[0], p2);
    retval[1] = modulo(p1[1], p2);
    retval[2] = modulo(p1[2], p2);
    return retval;
}

inline Position3 abs(Position3 const& v)
{
    T_ retval;
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
    T_ retval;
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

#define CREATE_VECTOR_INNER_TPL(__z__, __n__, __d__) \
    __d__[__n__] = BOOST_PP_CAT(p, __n__);

#define CREATE_VECTOR_TPL(__z__, __n__, __d__) \

template<typename T_> \
inline T_ create_vector( \
    BOOST_PP_ENUM_PARAMS(__n__, typename element_type_of<T_>::type const& p), \
    typename boost::enable_if<is_vector<T_, __n__> >::type* = 0) \
{ \
    T_ retval; \
    BOOST_PP_REPEAT_ ## __z__(__n__, CREATE_VECTOR_INNER_TPL, retval) \
    return retval; \
}

BOOST_PP_REPEAT_FROM_TO(2, CREATE_VECTOR_LIMIT_REPEAT, CREATE_VECTOR_TPL, )

#undef CREATE_VECTOR_TPL
#undef CREATE_VECTOR_INNER_TPL

template<typename T_>
inline bool is_cartesian_vector(Position3 const& vector)
{
    return (vector == create_vector<T_>(1, 0, 0) ||
            vector == create_vector<T_>(0, 1, 0) ||
            vector == create_vector<T_>(0, 0, 1) ||
            vector == create_vector<T_>(-1, 0, 0) ||
            vector == create_vector<T_>(0, -1, 0) ||
            vector == create_vector<T_>(0, 0, -1));
}

} // ecell4

#endif /* __LINEAR_ALGEBRA_HPP */
