#ifndef POSITION_HPP
#define POSITION_HPP

#include <ostream>
#include <functional>
#include <algorithm>
#include <cmath>
#include <boost/array.hpp>
#include <boost/range/begin.hpp>
#include <boost/range/end.hpp>

template<typename T_>
struct position: public boost::array<T_, 3>
{
    typedef boost::array<T_, 3> base_type;
    typedef typename base_type::value_type value_type;
    typedef typename base_type::size_type size_type;

    position()
    {
        (*this)[0] = 0;
        (*this)[1] = 0;
        (*this)[2] = 0;
    }

    position(const T_ (&a)[3]): base_type(
            *reinterpret_cast<const base_type*>(&a)) {}

    position(const T_ a[3]): base_type(
            *reinterpret_cast<const base_type*>(a)) {}

    position(const base_type& a): base_type(a) {}

    position(value_type x, value_type y, value_type z)
    {
        (*this)[0] = x;
        (*this)[1] = y;
        (*this)[2] = z;
    }

    value_type& x()
    {
        return (*this)[0];
    }

    const value_type& x() const
    {
        return (*this)[0];
    }

    value_type& y()
    {
        return (*this)[1];
    }

    const value_type& y() const
    {
        return (*this)[1];
    }

    value_type& z()
    {
        return (*this)[2];
    }

    const value_type& z() const
    {
        return (*this)[2];
    }


    const value_type distance_sq(const position& that) const
    {
        return std::pow((*this)[0] - that[0], 2)
            + std::pow((*this)[1] - that[1], 2)
            + std::pow((*this)[2] - that[2], 2);
    }

    const value_type distance(const position& that) const
    {
        return std::sqrt(distance_sq(that));
    }

    position operator+(const position& that) const
    {
        position retval;
        std::transform(
            boost::const_begin(*this), boost::const_end(*this),
            boost::const_begin(thabt), boost::begin(retval),
            std::plus<value_type>());
        return retval;
    }
};

template<typename Tstrm_, typename T_>
inline std::basic_ostream<Tstrm_>& operator<<(std::basic_ostream<Tstrm_>& strm,
        const position<T_>& v)
{
    strm << "(" << v.x() <<  ", " << v.y() <<  ", " << v.z() << ")";
    return strm;
}

#endif /* POSITION_HPP */
