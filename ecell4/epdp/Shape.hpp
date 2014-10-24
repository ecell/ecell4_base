#ifndef SHAPE_HPP
#define SHAPE_HPP

#include <boost/type_traits/remove_cv.hpp>
#include "Vector3.hpp"
#include "Position3Type.hpp"
#include "geometry.hpp"

template<typename Tshape_>
struct is_shape: public boost::mpl::false_ {};

template<typename Tobj_>
inline typename Tobj_::shape_type const& shape(Tobj_ const& obj)
{
    return obj.shape();
}

//XXX: See Shell.hpp and ParticleTraits.hpp
// template<typename Tobj_>
// inline typename Tobj_::shape_type& shape(Tobj_& obj)
// {
//     return obj.shape();
// }
//
// template<typename Tshape_>
// inline Tshape_ offset(Tshape_ const& shape, typename Tshape_::position_type off,
//         typename boost::enable_if<is_shape<Tshape_> >::type* = 0)
// {
//     Tshape_ retval(shape);
//     retval.position() += off;
//     return retval;
// }

template<typename T_>
struct shape_position_type
{
    struct argument_is_not_a_shape;
    static const std::size_t x = sizeof(argument_is_not_a_shape);
};

template<typename T_>
struct shape_length_type
{
    typedef typename element_type_of<typename shape_position_type<typename boost::remove_cv<T_>::type >::type>::type type;
};

template<typename T_>
inline typename shape_position_type<T_>::type const& shape_position(T_ const& shape,
        typename boost::enable_if<is_shape<T_> >::type* = 0)
{
    return shape.position();
}

template<typename T_>
inline typename shape_position_type<T_>::type& shape_position(T_& shape,
        typename boost::enable_if<is_shape<T_> >::type* = 0)
{
    return shape.position();
}

template<typename T_>
inline typename shape_length_type<T_>::type const& shape_size(T_ const& shape,
        typename boost::enable_if<is_shape<T_> >::type* = 0)
{
}

template<typename T_>
inline typename shape_length_type<T_>::type& shape_size(T_& shape,
        typename boost::enable_if<is_shape<T_> >::type* = 0)
{
}

template< typename T1_, typename T2_ >
inline typename shape_length_type<T1_>::type
distance_cyclic(
        T1_ const& p1, T2_ const& p2,
        typename shape_length_type<T1_>::type const& world_size,
        typename boost::enable_if<is_shape<T1_> >::type* = 0)
{
    return distance(p1, cyclic_transpose(p2, shape_position(p1), world_size));
}

template< typename T1_, typename T2_ >
inline typename shape_length_type<T1_>::type
distance_cyclic(
        T1_ const& p1, T2_ const& p2,
        typename shape_position_type<T1_>::type const& edge_lengths,
        typename boost::enable_if<is_shape<T1_> >::type* = 0)
{
    return distance(p1, cyclic_transpose(p2, shape_position(p1), edge_lengths));
}

template<typename T, typename Trng>
inline typename shape_position_type<T>::type
random_position(T const& shape, Trng const& rng)
{
    return random_position(shape, const_cast<Trng&>(rng));
}

#endif /* SHAPE_HPP */
