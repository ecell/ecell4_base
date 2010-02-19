#ifndef VECTOR3_HPP
#define VECTOR3_HPP

#include <ostream>
#include <functional>
#include <algorithm>
#include <cmath>
#include <gsl/gsl_pow_int.h>
#include <boost/mpl/bool.hpp>
#include <boost/type_traits/is_arithmetic.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/array.hpp>
#include <boost/range/begin.hpp>
#include <boost/range/end.hpp>

#if defined(HAVE_TR1_FUNCTIONAL)
#include <tr1/functional>
#elif defined(HAVE_STD_HASH)
#include <functional>
#elif defined(HAVE_BOOST_FUNCTIONAL_HASH_HPP)
#include <boost/functional/hash.hpp>
#endif

#include "utils/array_traits.hpp"

template<typename T_>
inline T_ add( T_ const& p1, T_ const& p2, typename boost::enable_if<boost::is_arithmetic<T_> >::type* = 0)
{
    return p1 + p2;
}

template<typename T_>
inline T_ subtract( T_ const& p1, T_ const& p2, typename boost::enable_if<boost::is_arithmetic<T_> >::type* = 0)
{
    return p1 - p2;
}

template<typename T_>
inline T_ multiply( T_ const& p1, T_ const& p2, typename boost::enable_if<boost::is_arithmetic<T_> >::type* = 0)
{
    return p1 * p2;
}

template<typename T_>
inline T_ divide( T_ const& p1, T_ const& p2, typename boost::enable_if<boost::is_arithmetic<T_> >::type* = 0)
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
inline T_ negate(T_ const& v, typename boost::enable_if<boost::is_arithmetic<T_> >::type* = 0)
{
    return -v;
}

template<typename T_>
struct is_vector3: public boost::mpl::false_ {};

template< typename T_ >
inline T_ add( T_ const& p1, T_ const& p2, typename boost::enable_if<is_vector3<T_> >::type* = 0 )
{
    T_ retval;
    retval[0] = add(p1[0], p2[0]);
    retval[1] = add(p1[1], p2[1]);
    retval[2] = add(p1[2], p2[2]);
    return retval;
}

template< typename T_ >
inline T_ subtract( T_ const& p1, T_ const& p2, typename boost::enable_if<is_vector3<T_> >::type* = 0 )
{
    T_ retval;
    retval[0] = subtract(p1[0], p2[0]);
    retval[1] = subtract(p1[1], p2[1]);
    retval[2] = subtract(p1[2], p2[2]);
    return retval;
}

template<typename T_>
inline T_ divide( T_ const& p1, typename T_::value_type p2, typename boost::enable_if<is_vector3<T_> >::type* = 0 )
{
    T_ retval;
    retval[0] = divide(p1[0], p2);
    retval[1] = divide(p1[1], p2);
    retval[2] = divide(p1[2], p2);
    return retval;
}

template<typename T_>
inline T_ multiply( T_ const& p1, typename T_::value_type p2, typename boost::enable_if<is_vector3<T_> >::type* = 0 )
{
    T_ retval;
    retval[0] = multiply(p1[0], p2);
    retval[1] = multiply(p1[1], p2);
    retval[2] = multiply(p1[2], p2);
    return retval;
}

template<typename T_>
inline T_ modulo( T_ const& p1, typename T_::value_type p2, typename boost::enable_if<is_vector3<T_> >::type* = 0 )
{
    T_ retval;
    retval[0] = modulo(p1[0], p2);
    retval[1] = modulo(p1[1], p2);
    retval[2] = modulo(p1[2], p2);
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
inline typename T_::value_type dot_product(T_ const& p1, T_ const& p2, typename boost::enable_if<is_vector3<T_> >::type* = 0)
{
    return multiply(p1[0], p2[0])
           + multiply(p1[1], p2[1])
           + multiply(p1[2], p2[2]);
}

template< typename T_ >
inline typename element_type_of< T_ >::type length_sq( T_ const& r )
{
    return gsl_pow_2( r[0] ) + gsl_pow_2( r[1] ) + gsl_pow_2( r[2] );
}

template< typename T_ >
inline typename element_type_of< T_ >::type length( T_ const& r )
{
    return std::sqrt( length_sq( r ) );
}

template<typename T_>
struct Vector3: public boost::array<T_, 3>
{
    typedef boost::array<T_, 3> base_type;
    typedef typename base_type::value_type value_type;
    typedef typename base_type::size_type size_type;

    Vector3& operator+=(Vector3 const& rhs)
    {
        *this = add(*this, rhs);
        return *this;
    }

    Vector3& operator-=(Vector3 const& rhs)
    {
        *this = subtract(*this, rhs);
        return *this;
    }

    template<typename TT_>
    Vector3& operator*=(TT_ const& rhs)
    {
        *this = multiply(*this, rhs);
        return *this;
    }

    template<typename TT_>
    Vector3& operator/=(TT_ const& rhs)
    {
        *this = divide(*this, rhs);
        return *this;
    }

    Vector3()
    {
        (*this)[0] = 0;
        (*this)[1] = 0;
        (*this)[2] = 0;
    }

    Vector3(const T_ (&a)[3]): base_type(
            *reinterpret_cast<const base_type*>(&a)) {}

    Vector3(const T_ a[3]): base_type(
            *reinterpret_cast<const base_type*>(a)) {}

    Vector3(const base_type& a): base_type(a) {}

    Vector3(value_type p0, value_type p1, value_type p2)
    {
        (*this)[0] = p0;
        (*this)[1] = p1;
        (*this)[2] = p2;
    }
};

template< typename T_ >
inline Vector3< T_ > operator+(Vector3< T_ > const& lhs, Vector3< T_ > const& rhs)
{
    return add( lhs, rhs );
}

template< typename T_ >
inline Vector3< T_ > operator-(Vector3< T_ > const& lhs, Vector3< T_ > const& rhs)
{
    return subtract( lhs, rhs );
}

template<typename T_>
inline Vector3<T_> operator/(Vector3<T_> const& lhs, T_ const& rhs)
{
    return divide(lhs, rhs);
}

template<typename T_>
inline Vector3<T_> operator*(Vector3<T_> const& lhs, T_ const& rhs)
{
    return multiply(lhs, rhs);
}

template<typename Tstrm_, typename Ttraits_, typename T_>
inline std::basic_ostream<Tstrm_, Ttraits_>&
operator<<(std::basic_ostream<Tstrm_, Ttraits_>& strm, const Vector3<T_>& v)
{
    strm << "(" << v[0] <<  ", " << v[1] <<  ", " << v[2] << ")";
    return strm;
}

template<typename T_>
struct is_vector3<Vector3<T_> >: public boost::mpl::true_ {};

template< typename T_ >
struct element_type_of< Vector3< T_ > >
{
    typedef T_ type;
};

#if defined(HAVE_TR1_FUNCTIONAL)
namespace std { namespace tr1 {
#elif defined(HAVE_STD_HASH)
namespace std {
#elif defined(HAVE_BOOST_FUNCTIONAL_HASH_HPP)
namespace boost {
#endif

template<typename T_>
struct hash<Vector3<T_> >
{
    typedef Vector3<T_> argument_type;

    std::size_t operator()(argument_type const& val)
    {
        return hash<typename argument_type::value_type>()(val[0]) ^
            hash<typename argument_type::value_type>()(val[1]) ^
            hash<typename argument_type::value_type>()(val[2]);
    }
};

#if defined(HAVE_TR1_FUNCTIONAL)
} } // namespace std::tr1
#elif defined(HAVE_STD_HASH)
} // namespace std
#elif defined(HAVE_BOOST_FUNCTIONAL_HASH_HPP)
} // namespace boost
#endif

#endif /* VECTOR3_HPP */
