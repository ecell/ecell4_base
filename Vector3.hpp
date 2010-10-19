#ifndef VECTOR3_HPP
#define VECTOR3_HPP

#include <ostream>
#include <iomanip>
#include <functional>
#include <algorithm>

#if defined(HAVE_TR1_FUNCTIONAL)
#include <tr1/functional>
#elif defined(HAVE_STD_HASH)
#include <functional>
#elif defined(HAVE_BOOST_FUNCTIONAL_HASH_HPP)
#include <boost/functional/hash.hpp>
#endif

#include <boost/array.hpp>
#include "utils/array_traits.hpp"
#include "linear_algebra.hpp"

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
    strm << std::setprecision(12) << "(" << v[0] <<  ", " << v[1] <<  ", " << v[2] << ")";
    return strm;
}

template<typename T_>
struct is_vector<Vector3<T_>, 3>: public boost::mpl::true_ {};

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
