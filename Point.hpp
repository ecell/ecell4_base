#ifndef POINT_HPP
#define POINT_HPP

#include "Shape.hpp"
#include "Vector3.hpp"

template<typename T_>
struct shape_position_type<Vector3<T_> >
{
    typedef Vector3<T_> type;
};

template<typename T_>
struct shape_length_type<Vector3<T_> >
{
    typedef T_ type;
};

template<typename T_>
inline Vector3<T_> shape_position(Vector3<T_> const& v)
{
    return v;
}

#endif /* POINT_HPP */
