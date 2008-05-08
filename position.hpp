#ifndef POSITION_HPP
#define POSITION_HPP

#include <ostream>
#include <boost/array.hpp>

template<typename T_>
struct position: public boost::array<T_, 3>
{
    typedef boost::array<T_, 3> base_type;

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

    position(T_ x, T_ y, T_ z)
    {
        (*this)[0] = x;
        (*this)[1] = y;
        (*this)[2] = z;
    }

    T_& x() {
        return (*this)[0];
    }

    const T_& x() const {
        return (*this)[0];
    }

    T_& y() {
        return (*this)[1];
    }

    const T_& y() const {
        return (*this)[1];
    }

    T_& z() {
        return (*this)[2];
    }

    const T_& z() const {
        return (*this)[2];
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
