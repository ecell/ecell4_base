#ifndef ECELL4_MAKE_VISITOR
#define ECELL4_MAKE_VISITOR
#include <boost/static_assert.hpp>
#include <boost/preprocessor.hpp>
#include <boost/function.hpp>
#include <boost/bind.hpp>

#ifndef ECELL4_MAKE_VISITOR_MAX_SIZE
#define ECELL4_MAKE_VISITOR_MAX_SIZE 10
#endif//ECELL4_MAKE_VISITOR_MAX_SIZE

#define ECELL4_MAKE_VISITOR_MAX_INDEX BOOST_PP_ADD(ECELL4_MAKE_VISITOR_MAX_SIZE, 1)

namespace ecell4
{

//TODO:
// - assert if not all the result_type is same

template<typename T>
struct function_traits
{
    typedef typename T::result_type result_type;
};

template<typename R, typename T>
struct function_traits<boost::function<R(T)> >
{
    typedef R result_type;
};

namespace detail
{

// ---------------------------- aggregate_functions ----------------------------

#define EXPAND_BASECLASS_INITIALIZER(z, N, NAME)\
    BOOST_PP_CAT(NAME, N)(BOOST_PP_CAT(t, N))\
    /**/

#define EXPAND_USING_OPERATORS(z, N, NAME)\
    using BOOST_PP_CAT(NAME, N)::operator();\
    /**/

#define EXPAND_AGGREGATE_FUNCTIONS(z, N, DUMMY)\
template<BOOST_PP_ENUM_PARAMS(N, typename T)>\
struct BOOST_PP_CAT(aggregate_functions, N)\
    :BOOST_PP_ENUM_PARAMS(N, T)\
{\
    typedef typename function_traits<T0>::result_type result_type;\
\
    BOOST_PP_CAT(aggregate_functions, N)(BOOST_PP_ENUM_BINARY_PARAMS(N, T, t))\
    :BOOST_PP_ENUM(N, EXPAND_BASECLASS_INITIALIZER, T)\
    {}\
\
    BOOST_PP_REPEAT(N, EXPAND_USING_OPERATORS, T)\
};\
/**/

BOOST_PP_REPEAT_FROM_TO(1, ECELL4_MAKE_VISITOR_MAX_INDEX, EXPAND_AGGREGATE_FUNCTIONS, dummy)

#undef EXPAND_CONSTRUCTOR_ARGUMENTS
#undef EXPAND_BASECLASS_INITIALIZER
#undef EXPAND_AGGREGATE_FUNCTIONS

} // detail

#define EXPAND_MAKE_VISITOR(z, N, DUMMY)\
template<BOOST_PP_ENUM_PARAMS(N, typename T)>\
inline BOOST_PP_CAT(detail::aggregate_functions, N)<BOOST_PP_ENUM_PARAMS(N, T)>\
make_visitor(BOOST_PP_ENUM_BINARY_PARAMS(N, T, t))\
{\
    return BOOST_PP_CAT(detail::aggregate_functions, N)<BOOST_PP_ENUM_PARAMS(N, T)>(\
            BOOST_PP_ENUM_PARAMS(N, t)\
            );\
}\
/**/

BOOST_PP_REPEAT_FROM_TO(1, ECELL4_MAKE_VISITOR_MAX_INDEX, EXPAND_MAKE_VISITOR, dummy)

#undef EXPAND_MAKE_VISITOR

template<typename R, typename T>
inline boost::function<R(T)>
resolve(R(*fptr)(T))
{
    return boost::function<R(T)>(fptr);
}

template<typename R, typename T, class C>
inline boost::function<R(T)>
resolve(R(C::*fptr)(T), C* cptr)
{
    return boost::function<R(T)>(boost::bind(fptr, cptr, _1));
}

} // ecell4

#undef ECELL4_MAKE_VISITOR_MAX_INDEX
#endif// ECELL4_MAKE_OVERLOAD
