#ifndef LINEAR_ALGEBRA_HPP
#define LINEAR_ALGEBRA_HPP

#include <algorithm>
#include <cmath>
#include <gsl/gsl_pow_int.h>
#include <boost/preprocessor/tuple/elem.hpp>
#include <boost/preprocessor/repetition/repeat.hpp>
#include <boost/preprocessor/repetition/enum_params.hpp>
#include <boost/preprocessor/repetition/enum_trailing_params.hpp>
#include <boost/mpl/bool.hpp>
#include <boost/mpl/and.hpp>
#include <boost/type_traits/is_arithmetic.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/array.hpp>
#include <boost/multi_array.hpp>
#include <boost/preprocessor/repetition/repeat_from_to.hpp>

#include "utils/array_traits.hpp"

#include <ecell4/core/functions.hpp>
using ecell4::pow_2;

#define CREATE_VECTOR_LIMIT_REPEAT 16
#define POPULATE_MATRIX_BY_VECTORS_LIMIT_REPEAT 16

template<typename T_, std::size_t N_>
struct is_vector: public boost::mpl::false_ {};

template<typename T_, std::size_t N_>
struct is_vector<boost::array<T_, N_>, N_>: public boost::mpl::true_ {};

template<typename T_, std::size_t N_>
struct is_matrix: public boost::mpl::false_ {};

template<typename T_, std::size_t N_, typename Talloc_>
struct is_matrix<boost::multi_array<T_, N_, Talloc_>, N_>: public boost::mpl::true_ {};

template<typename T_, std::size_t N_>
struct is_matrix<boost::detail::multi_array::sub_array<T_, N_>, N_>: public boost::mpl::true_ {};

template<typename T_, std::size_t N_>
struct is_matrix<boost::detail::multi_array::const_sub_array<T_, N_>, N_>: public boost::mpl::true_ {};

template<typename T_, std::size_t N_>
struct is_matrix<boost::detail::multi_array::multi_array_view<T_, N_>, N_>: public boost::mpl::true_ {};

template<typename T_, std::size_t N_>
struct is_matrix<boost::detail::multi_array::const_multi_array_view<T_, N_>, N_>: public boost::mpl::true_ {};

template<typename T_, std::size_t N_>
struct is_matrix<boost::multi_array_ref<T_, N_>, N_>: public boost::mpl::true_ {};

template<typename T_, std::size_t N_, typename Tptr_>
struct is_matrix<boost::const_multi_array_ref<T_, N_, Tptr_>, N_>: public boost::mpl::true_ {};

template<typename T_, std::size_t N1_, std::size_t N2_>
struct is_matrix<boost::array<boost::array<T_, N1_>, N2_>, 2>: public boost::mpl::true_ {};

template<typename T_, std::size_t N1_, std::size_t N2_>
struct is_matrix<T_[N1_][N2_], 2>: public boost::mpl::true_ {};

template<typename T_>
struct is_scalar: public boost::is_arithmetic<T_> {};

template<typename T_>
struct is_vector2: public is_vector<T_, 2> {};

template<typename T_>
struct is_vector3: public is_vector<T_, 3> {};

template<typename T_>
struct matrix_adapter
{
};

template<typename T_, std::size_t N_, typename Talloc_>
struct matrix_adapter<boost::multi_array<T_, N_, Talloc_> >
{
    typedef boost::multi_array<T_, N_, Talloc_> matrix_type;
    typedef typename matrix_type::size_type matrix_size_type;

    static matrix_size_type get_extent(matrix_type const& first,
                                       std::size_t idx)
    {
        return first.shape()[idx];
    }
};

template<typename T_, std::size_t N_>
struct matrix_adapter<boost::detail::multi_array::sub_array<T_, N_> >
{
    typedef boost::detail::multi_array::sub_array<T_, N_> matrix_type;
    typedef typename matrix_type::size_type matrix_size_type;

    static matrix_size_type get_extent(matrix_type const& first,
                                       std::size_t idx)
    {
        return first.shape()[idx];
    }
};

template<typename T_, std::size_t N_>
struct matrix_adapter<boost::detail::multi_array::const_sub_array<T_, N_> >
{
    typedef boost::detail::multi_array::const_sub_array<T_, N_> matrix_type;
    typedef typename matrix_type::size_type matrix_size_type;

    static matrix_size_type get_extent(matrix_type const& first,
                                       std::size_t idx)
    {
        return first.shape()[idx];
    }
};

template<typename T_, std::size_t N_>
struct matrix_adapter<boost::detail::multi_array::multi_array_view<T_, N_> >
{
    typedef boost::detail::multi_array::multi_array_view<T_, N_> matrix_type;
    typedef typename matrix_type::size_type matrix_size_type;

    static matrix_size_type get_extent(matrix_type const& first,
                                       std::size_t idx)
    {
        return first.shape()[idx];
    }
};

template<typename T_, std::size_t N_>
struct matrix_adapter<boost::detail::multi_array::const_multi_array_view<T_, N_> >
{
    typedef boost::detail::multi_array::const_multi_array_view<T_, N_> matrix_type;
    typedef typename matrix_type::size_type matrix_size_type;

    matrix_size_type get_extent(matrix_type const& first,
                                std::size_t idx)
    {
        return first.shape()[idx];
    }
};

template<typename T_, std::size_t N_>
struct matrix_adapter<boost::multi_array_ref<T_, N_> >
{
    typedef boost::multi_array_ref<T_, N_> matrix_type;
    typedef typename matrix_type::size_type matrix_size_type;

    static matrix_size_type get_extent(matrix_type const& first,
                                       std::size_t idx)
    {
        return first.shape()[idx];
    }
};

template<typename T_, std::size_t N_, typename Tptr_>
struct matrix_adapter<boost::const_multi_array_ref<T_, N_, Tptr_> >
{
    typedef boost::const_multi_array_ref<T_, N_, Tptr_> matrix_type;
    typedef typename matrix_type::size_type matrix_size_type;

    static matrix_size_type get_extent(matrix_type const& first,
                                       std::size_t idx)
    {
        return first.shape()[idx];
    }
};

template<typename T_, std::size_t N1_, std::size_t N2_>
struct matrix_adapter<T_[N1_][N2_]>
{
    typedef T_ matrix_type[N1_][N2_];
    typedef std::size_t matrix_size_type;

    static matrix_size_type get_extent(matrix_type const& first,
                                       std::size_t idx)
    {
        switch (idx)
        {
        case 0:
            return N1_;
        case 1:
            return N2_;
        default:
            throw std::out_of_range("index out of range");
        }
    }
};

template<typename Tmat>
inline std::size_t matrix_extent(Tmat const& mat, std::size_t dim)
{
    return matrix_adapter<Tmat>::get_extent(mat, dim);
}

template<typename T_>
inline T_ add( T_ const& p1, T_ const& p2, typename boost::enable_if<is_scalar<T_> >::type* = 0)
{
    return p1 + p2;
}

template<typename T_>
inline T_ subtract( T_ const& p1, T_ const& p2, typename boost::enable_if<is_scalar<T_> >::type* = 0)
{
    return p1 - p2;
}

template<typename T_>
inline T_ multiply( T_ const& p1, T_ const& p2, typename boost::enable_if<is_scalar<T_> >::type* = 0)
{
    return p1 * p2;
}

template<typename T_>
inline T_ divide( T_ const& p1, T_ const& p2, typename boost::enable_if<is_scalar<T_> >::type* = 0)
{
    return p1 / p2;
}

template<typename T_>
inline T_ modulo( T_ const& p1, T_ const& p2 )
{
    T_ r = p1 % p2;
    if (r != 0 && (r > 0) == (p2 < 0))
        r += p2;
    return r;
}

template<>
inline float modulo( float const& p1, float const& p2 )
{
    float r = std::fmod(p1, p2);
    if (r != 0 && (r > 0) == (p2 < 0))
        r += p2;
    return r;
}

template<>
inline double modulo( double const& p1, double const& p2 )
{
    double r = std::fmod(p1, p2);
    if (r != 0 && (r > 0) == (p2 < 0))
        r += p2;
    return r;
}

template<typename T_>
inline T_ negate(T_ const& v, typename boost::enable_if<is_scalar<T_> >::type* = 0)
{
    return -v;
}

template<typename T_>
inline T_ abs(T_ const& v, typename boost::enable_if<is_scalar<T_> >::type* = 0)
{
    return std::fabs(v);
}

template<typename T1_, typename T2_>
inline T1_ add(T1_ const& p1, T2_ const& p2, typename boost::enable_if<boost::mpl::and_<is_vector3<T1_>, is_vector3<T2_> > >::type* = 0)
{
    T1_ retval;
    retval[0] = add(p1[0], p2[0]);
    retval[1] = add(p1[1], p2[1]);
    retval[2] = add(p1[2], p2[2]);
    return retval;
}

template<typename T1_, typename T2_>
inline T1_ subtract(T1_ const& p1, T2_ const& p2, typename boost::enable_if<boost::mpl::and_<is_vector3<T1_>, is_vector3<T2_> > >::type* = 0)
{
    T1_ retval;
    retval[0] = subtract(p1[0], p2[0]);
    retval[1] = subtract(p1[1], p2[1]);
    retval[2] = subtract(p1[2], p2[2]);
    return retval;
}

template<typename T_>
inline T_ divide(T_ const& p1, typename element_type_of<T_>::type const& p2, typename boost::enable_if<is_vector3<T_> >::type* = 0)
{
    T_ retval;
    retval[0] = divide(p1[0], p2);
    retval[1] = divide(p1[1], p2);
    retval[2] = divide(p1[2], p2);
    return retval;
}

template<typename T_>
inline T_ multiply(T_ const& p1, typename element_type_of<T_>::type const& p2, typename boost::enable_if<is_vector3<T_> >::type* = 0)
{
    T_ retval;
    retval[0] = multiply(p1[0], p2);
    retval[1] = multiply(p1[1], p2);
    retval[2] = multiply(p1[2], p2);
    return retval;
}

template<typename T_, typename M_>
inline T_ multiply(T_ const& p1, M_ const& p2, typename boost::enable_if<
    boost::mpl::and_<is_vector3<T_>, is_matrix<M_, 2> > >::type* = 0)
{
    BOOST_ASSERT(matrix_extent(p2, 0) == 3 && matrix_extent(p2, 1) == 3);
    T_ retval;
    retval[0] = multiply(p1[0], p2[0][0])
                + multiply(p1[1], p2[1][0])
                + multiply(p1[2], p2[2][0]);
    retval[1] = multiply(p1[0], p2[0][1])
                + multiply(p1[1], p2[1][1])
                + multiply(p1[2], p2[2][1]);
    retval[2] = multiply(p1[0], p2[0][2])
                + multiply(p1[1], p2[1][2])
                + multiply(p1[2], p2[2][2]);
    return retval;
}

template<typename M_, typename T_>
inline T_ multiply(M_ const& p1, T_ const& p2, typename boost::enable_if<
    boost::mpl::and_<is_matrix<M_, 2>, is_vector3<T_> > >::type* = 0)
{
    BOOST_ASSERT(matrix_extent(p1, 0) == 3 && matrix_extent(p1, 1) == 3);
    T_ retval;
    retval[0] = multiply(p1[0][0], p2[0])
                + multiply(p1[0][1], p2[1])
                + multiply(p1[0][2], p2[2]);
    retval[1] = multiply(p1[1][0], p2[0])
                + multiply(p1[1][1], p2[1])
                + multiply(p1[1][2], p2[2]);
    retval[2] = multiply(p1[2][0], p2[0])
                + multiply(p1[2][1], p2[1])
                + multiply(p1[2][2], p2[2]);
    return retval;
}

template<typename T_>
inline T_ modulo(T_ const& p1, typename element_type_of<T_>::type const& p2, typename boost::enable_if<is_vector3<T_> >::type* = 0)
{
    T_ retval;
    retval[0] = modulo(p1[0], p2);
    retval[1] = modulo(p1[1], p2);
    retval[2] = modulo(p1[2], p2);
    return retval;
}

template<typename T_>
inline T_ modulo(T_ const& p1, T_ const& p2, typename boost::enable_if<is_vector3<T_> >::type* = 0)
{
    T_ retval;
    retval[0] = modulo(p1[0], p2[0]);
    retval[1] = modulo(p1[1], p2[1]);
    retval[2] = modulo(p1[2], p2[2]);
    return retval;
}

template<typename T_>
inline T_ negate(T_ const& v, typename boost::enable_if<is_vector3<T_> >::type* = 0)
{
    T_ retval;
    retval[0] = -v[0];
    retval[1] = -v[1];
    retval[2] = -v[2];
    return retval;
}

template<typename T_>
inline T_ abs(T_ const& v, typename boost::enable_if<is_vector3<T_> >::type* = 0)
{
    T_ retval;
    retval[0] = abs(v[0]);
    retval[1] = abs(v[1]);
    retval[2] = abs(v[2]);
    return retval;
}

template<typename T_>
inline typename element_type_of<T_>::type dot_product(T_ const& p1, T_ const& p2, typename boost::enable_if<is_vector3<T_> >::type* = 0)
{
    return multiply(p1[0], p2[0])
           + multiply(p1[1], p2[1])
           + multiply(p1[2], p2[2]);
}

template<typename T_>
inline T_ cross_product(T_ const& p1, T_ const& p2, typename boost::enable_if<is_vector3<T_> >::type* = 0)
{
    T_ retval;
    retval[0] = subtract(multiply(p1[1], p2[2]), multiply(p1[2], p2[1]));
    retval[1] = subtract(multiply(p1[2], p2[0]), multiply(p1[0], p2[2]));
    retval[2] = subtract(multiply(p1[0], p2[1]), multiply(p1[1], p2[0]));
    return retval;
}

template<typename T_>
inline typename element_type_of<T_>::type length_sq(T_ const& r, typename boost::enable_if<is_vector2<T_> >::type* = 0)
{
    return pow_2(r[0]) + pow_2(r[1]);
}

template< typename T_ >
inline typename element_type_of< T_ >::type length_sq(T_ const& r, typename boost::enable_if<is_vector3<T_> >::type* = 0)
{
    return pow_2(r[0]) + pow_2(r[1]) + pow_2(r[2]);
}

template< typename T_ >
inline typename element_type_of< T_ >::type length(T_ const& r)
{
    return std::sqrt(length_sq(r));
}

#define CREATE_VECTOR_INNER_TPL(__z__, __n__, __d__) \
    __d__[__n__] = BOOST_PP_CAT(p, __n__);

#define CREATE_VECTOR_TPL(__z__, __n__, __d__) \
template<typename T_> \
inline T_ create_vector(\
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
inline bool is_cartesian_versor(T_ const& vector, typename boost::enable_if<is_vector3<T_> >::type* = 0)
{
    return (vector == create_vector<T_>(1, 0, 0) ||
            vector == create_vector<T_>(0, 1, 0) ||
            vector == create_vector<T_>(0, 0, 1) ||
            vector == create_vector<T_>(-1, 0, 0) ||
            vector == create_vector<T_>(0, -1, 0) ||
            vector == create_vector<T_>(0, 0, -1));
}

#endif /* LINEAR_ALGEBRA_HPP */
