#ifndef ECELL4_MAKE_VISITOR_HPP
#define ECELL4_MAKE_VISITOR_HPP
#include <boost/preprocessor.hpp>
#include <boost/typeof/typeof.hpp>
#include <boost/type_traits.hpp>
#include <boost/mpl/assert.hpp>
#include <boost/mpl/logical.hpp>
#include <boost/mpl/and.hpp>
#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <boost/variant.hpp>

#ifndef ECELL4_MAKE_VISITOR_MAX_SIZE
#define ECELL4_MAKE_VISITOR_MAX_SIZE BOOST_VARIANT_LIMIT_TYPES
#endif//ECELL4_MAKE_VISITOR_MAX_SIZE

#define ECELL4_MAKE_VISITOR_MAX_INDEX BOOST_PP_ADD(ECELL4_MAKE_VISITOR_MAX_SIZE, 1)

namespace ecell4
{

template<typename T>
struct result_type_of
{
    typedef typename T::result_type type;
};

namespace detail
{

#define ECELL4_AGGREGATE_FUNCTION_BASECLASS_INITIALIZER(z, N, NAME)\
    BOOST_PP_CAT(NAME, N)(BOOST_PP_CAT(t, N))\
    /**/

#define ECELL4_AGGREGATE_FUNCTION_USING_OPERATORS(z, N, NAME)\
    using BOOST_PP_CAT(NAME, N)::operator();\
    /**/

//XXX error: redeclaration of
//           `boost_static_assert_typedef_<line num @ the macro expanded>`
// #define ECELL4_CHECK_RESULT_TYPE_IS_SAME(z, N, data)\
//     BOOST_STATIC_ASSERT((boost::is_same<\
//         typename result_type_of<BOOST_PP_CAT(T, N)>::type, result_type>::type::value));
//     BOOST_PP_REPEAT_FROM_TO(1, N, ECELL4_CHECK_RESULT_TYPE_IS_SAME, DUMMY)\
// #undef ECELL4_CHECK_RESULT_TYPE_IS_SAME

#define ECELL4_AGGREGATE_FUNCTIONS(z, N, DUMMY)\
template<BOOST_PP_ENUM_PARAMS(N, typename T)>\
struct BOOST_PP_CAT(aggregate_functions, N)\
    :BOOST_PP_ENUM_PARAMS(N, T)\
{\
    typedef typename result_type_of<T0>::type result_type;\
\
    BOOST_PP_CAT(aggregate_functions, N)(BOOST_PP_ENUM_BINARY_PARAMS(N, T, t))\
    :BOOST_PP_ENUM(N, ECELL4_AGGREGATE_FUNCTION_BASECLASS_INITIALIZER, T)\
    {}\
\
    BOOST_PP_REPEAT(N, ECELL4_AGGREGATE_FUNCTION_USING_OPERATORS, T)\
};\
/**/

BOOST_PP_REPEAT_FROM_TO(1, ECELL4_MAKE_VISITOR_MAX_INDEX, ECELL4_AGGREGATE_FUNCTIONS, dummy)

#undef ECELL4_AGGREGATE_FUNCTION_BASECLASS_INITIALIZER
#undef ECELL4_AGGREGATE_FUNCTION_USING_OPERATORS
#undef ECELL4_AGGREGATE_FUNCTIONS

} // detail

#define ECELL4_MAKE_VISITOR(z, N, DUMMY)\
template<BOOST_PP_ENUM_PARAMS(N, typename T)>\
inline BOOST_PP_CAT(detail::aggregate_functions, N)<BOOST_PP_ENUM_PARAMS(N, T)>\
make_visitor(BOOST_PP_ENUM_BINARY_PARAMS(N, T, t))\
{\
    return BOOST_PP_CAT(detail::aggregate_functions, N)<BOOST_PP_ENUM_PARAMS(N, T)>(\
            BOOST_PP_ENUM_PARAMS(N, t)\
            );\
}\
/**/

BOOST_PP_REPEAT_FROM_TO(1, ECELL4_MAKE_VISITOR_MAX_INDEX, ECELL4_MAKE_VISITOR, dummy)

#undef ECELL4_MAKE_VISITOR_MAX_INDEX
#undef ECELL4_MAKE_VISITOR

template<typename R, typename T>
inline boost::visitor_ptr_t<T, R>
resolve(R(*fptr)(T))
{
    return boost::visitor_ptr(fptr);
}

//TODO: remove wrapper `boost::function` using some technique like
//      `decltype(boost::bind(declval<R(C::*)(T)>(), declval<C*>()))`
template<typename R, typename T, class C>
inline boost::function<R(T)>
resolve(R(C::*fptr)(T), C* cptr)
{
    return boost::function<R(T)>(boost::bind(fptr, cptr, _1));
}

} // ecell4
#endif// ECELL4_MAKE_VISITOR_HPP
