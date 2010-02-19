#ifndef SHAPE_HPP
#define SHAPE_HPP

#include "Vector3.hpp"

template<typename Tobj_>
inline typename Tobj_::shape_type const& shape(Tobj_ const& obj)
{
    return obj.shape();
}

template<typename Tobj_>
inline typename Tobj_::shape_type& shape(Tobj_& obj)
{
    return obj.shape();
}

template<typename Tshape_>
inline Tshape_ offset(Tshape_ const& shape, typename Tshape_::position_type off)
{
    Tshape_ retval(shape);
    retval.position() += off;
    return retval;
}

template<typename Tshape_>
struct is_shape: public boost::mpl::false_ {};

template<typename T_>
struct shape_position_type {};

template<typename T_>
struct shape_length_type
{
    typedef typename element_type_of<typename shape_position_type<T_>::type>::type type;
};

template<typename T_>
inline typename shape_position_type<T_>::type shape_position(T_ const& shape)
{
    return shape.position();
}

template< typename T1_, typename T2_ >
inline typename shape_length_type<T1_>::type
distance_cyclic(
        T1_ const& p1, T2_ const& p2,
        typename shape_length_type<T1_>::type const& world_size)
{
    return distance(p1, cyclic_transpose(p2, shape_position(p1), world_size));
}

template< typename T1_, typename T2_ >
inline typename shape_length_type<T1_>::type
distance_sq_cyclic(
        T1_ const& p1, T2_ const& p2,
        typename shape_length_type<T1_>::type const& world_size)
{
    return distance_sq(p1, cyclic_transpose(p2, shape_position(p1), world_size));
}


#endif /* SHAPE_HPP */
