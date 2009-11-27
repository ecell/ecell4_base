#ifndef FUN_WRAPPERS_HPP
#define FUN_WRAPPERS_HPP

#include <boost/utility/enable_if.hpp>
#include "utils/fun_composition.hpp"

template < typename T_ >
struct delete_ptr
{
    typedef void result_type;
    typedef T_* argument_type;

    void operator()( T_* ptr )
    {
        delete ptr;
    }
};

template<typename T_, typename Targ_>
struct reinterpret_caster
{
    T_ operator()(Targ_ const& v)
    {
        return reinterpret_cast<T_>(v);
    }
};

template<typename T_, typename Targ_>
struct reinterpret_caster<T_&, Targ_&>
{
    T_& operator()(Targ_& v)
    {
        return reinterpret_cast<T_&>(v);
    }
};

template<typename T_, typename Targ_>
inline T_ reinterpret_cast_wrapper(Targ_ const& v, typename boost::disable_if<boost::is_reference<T_> >::type* = 0)
{
    return reinterpret_caster<T_, Targ_>()(v);
}

template<typename T_, typename Targ_>
inline T_& reinterpret_cast_wrapper(Targ_& v)
{
    return reinterpret_caster<T_&, Targ_&>()(v);
}

template<typename T_, typename Targ_>
struct dynamic_caster
{
    T_ operator()(Targ_ const& v)
    {
        return dynamic_cast<T_>(v);
    }
};

template<typename T_, typename Targ_>
struct dynamic_caster<T_&, Targ_&>
{
    T_& operator()(Targ_& v)
    {
        return dynamic_cast<T_&>(v);
    }
};

template<typename T_, typename Targ_>
inline T_ dynamic_cast_wrapper(Targ_ const& v, typename boost::disable_if<boost::is_reference<T_> >::type* = 0)
{
    return dynamic_caster<T_, Targ_>()(v);
}

template<typename T_, typename Targ_>
inline T_& dynamic_cast_wrapper(Targ_& v)
{
    return dynamic_caster<T_&, Targ_&>()(v);
}

#endif /* FUN_WRAPPERS_HPP */
