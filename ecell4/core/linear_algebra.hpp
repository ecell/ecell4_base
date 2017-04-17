#ifndef ECELL4_LINEAR_ALGEBRA_HPP
#define ECELL4_LINEAR_ALGEBRA_HPP

#include <algorithm>

#include "Real3.hpp"

namespace ecell4
{

/*

#define CREATE_VECTOR_INNER_TPL(__z__, __n__, __d__) \
    __d__[__n__] = BOOST_PP_CAT(p, __n__);

#define CREATE_VECTOR_TPL(__z__, __n__, __d__) \

template<typename T_> \
inline T_ create_vector( \
    BOOST_PP_ENUM_PARAMS(__n__, const typename element_type_of<T_>::type& p), \
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
inline bool is_cartesian_vector(const Real3& vector)
{
    return (vector == create_vector<T_>(1, 0, 0) ||
            vector == create_vector<T_>(0, 1, 0) ||
            vector == create_vector<T_>(0, 0, 1) ||
            vector == create_vector<T_>(-1, 0, 0) ||
            vector == create_vector<T_>(0, -1, 0) ||
            vector == create_vector<T_>(0, 0, -1));
}

*/

} // ecell4

#endif /* ECELL4_LINEAR_ALGEBRA_HPP */
