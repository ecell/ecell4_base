#ifndef ARRAY_HELPER_HPP
#define ARRAY_HELPER_HPP

#include <boost/array.hpp>
#include <boost/preprocessor/config/limits.hpp>
#include <boost/preprocessor/repetition/repeat.hpp>
#include <boost/preprocessor/repetition/enum_params.hpp>
#include <boost/preprocessor/cat.hpp>
#include <boost/call_traits.hpp>

#define ARRAY_HELPER_INNER_TPL(__z__, __n__, __d__) \
    __d__[__n__] = BOOST_PP_CAT(p, __n__);

#define ARRAY_HELPER_TPL(__z__, __n__, __d__) \
template<typename T_> \
inline ::boost::array<T_, __n__> array_gen(\
        BOOST_PP_ENUM_PARAMS(__n__, T_ const& p)) \
{ \
    ::boost::array<T_, __n__> retval; \
    BOOST_PP_REPEAT_ ## __z__(__n__, ARRAY_HELPER_INNER_TPL, retval) \
    return retval; \
}

BOOST_PP_REPEAT(BOOST_PP_LIMIT_REPEAT, ARRAY_HELPER_TPL, )

#undef ARRAY_HELPER_TPL
#undef ARRAY_HELPER_INNER_TPL


#endif /* ARRAY_HELPER_HPP */
