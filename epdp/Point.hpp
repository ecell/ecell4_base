#ifndef POINT_HPP
#define POINT_HPP

#include "Shape.hpp"
#include "Vector3.hpp"

#include "Position3Type.hpp"

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


template<>
struct shape_position_type<ecell4::Position3>
{
    typedef ecell4::Position3 type;
};

template<>
struct shape_length_type<ecell4::Position3>
{
    typedef ecell4::Position3::value_type type;
};

inline ecell4::Position3 shape_position(ecell4::Position3 const &v)
{
    return v;
}
#endif /* POINT_HPP */
