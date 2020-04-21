#ifndef LINEAR_ALGEBRA_HPP
#define LINEAR_ALGEBRA_HPP

#include <algorithm>
#include <type_traits>
#include <cmath>
#include <gsl/gsl_pow_int.h>
#include <boost/mpl/and.hpp>
#include <array>

#include "utils/array_traits.hpp"

#include <ecell4/core/functions.hpp>

namespace ecell4
{
namespace egfrd
{

using ::ecell4::pow_2;

// ----------------------------------------------------------------------------
// is_vector

template<typename T_, std::size_t N_>
struct is_vector: public std::false_type {};

template<typename T_, std::size_t N_>
struct is_vector<std::array<T_, N_>, N_>: public std::true_type {};

// ----------------------------------------------------------------------------
// is_matrix

template<typename T_, std::size_t N_>
struct is_matrix: public std::false_type {};

template<typename T_, std::size_t N1_, std::size_t N2_>
struct is_matrix<T_[N1_][N2_], 2>: public std::true_type {};

// ----------------------------------------------------------------------------
// is_scalar

template<typename T_>
struct is_scalar: public boost::is_arithmetic<T_> {};

// ----------------------------------------------------------------------------
// is_vector3

template<typename T_>
struct is_vector3: public is_vector<T_, 3> {};

// ----------------------------------------------------------------------------
// matrix_adapter

template<typename T_>
struct matrix_adapter{};

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

// ----------------------------------------------------------------------------
// matrix_adapter

template<typename T_>
inline T_ add( T_ const& p1, T_ const& p2, typename std::enable_if<is_scalar<T_>::value>::type* = 0)
{
    return p1 + p2;
}

template<typename T_>
inline T_ subtract( T_ const& p1, T_ const& p2, typename std::enable_if<is_scalar<T_>::value>::type* = 0)
{
    return p1 - p2;
}

template<typename T_>
inline T_ multiply( T_ const& p1, T_ const& p2, typename std::enable_if<is_scalar<T_>::value>::type* = 0)
{
    return p1 * p2;
}

template<typename T_>
inline T_ divide( T_ const& p1, T_ const& p2, typename std::enable_if<is_scalar<T_>::value>::type* = 0)
{
    return p1 / p2;
}

template<typename T_>
inline T_ modulo( T_ const& p1, T_ const& p2 )
{
    T_ r = p1 % p2;
    if (r != 0 && (r > 0) == (p2 < 0))
    {
        r += p2;
    }
    return r;
}

template<>
inline float modulo( float const& p1, float const& p2 )
{
    float r = std::fmod(p1, p2);
    if (r != 0 && (r > 0) == (p2 < 0))
    {
        r += p2;
    }
    return r;
}

template<>
inline double modulo( double const& p1, double const& p2 )
{
    double r = std::fmod(p1, p2);
    if (r != 0 && (r > 0) == (p2 < 0))
    {
        r += p2;
    }
    return r;
}

template<typename T_, typename M_>
inline T_ multiply(T_ const& p1, M_ const& p2, typename std::enable_if<
    boost::mpl::and_<is_vector3<T_>, is_matrix<M_, 2> >::value>::type* = 0)
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
inline T_ multiply(M_ const& p1, T_ const& p2, typename std::enable_if<
    boost::mpl::and_<is_matrix<M_, 2>, is_vector3<T_> >::value>::type* = 0)
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

template<typename T, typename ... Args>
typename std::enable_if<is_vector<T, sizeof...(Args)>::value, T>::type
create_vector(Args&& ... args)
{
    return T(std::forward<Args>(args)...);
}

} // egfrd
} // ecell4
#endif /* LINEAR_ALGEBRA_HPP */
