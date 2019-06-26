#ifndef ECELL4_QUANTITY_HPP
#define ECELL4_QUANTITY_HPP

#include <ostream>
#include <string>
#include <boost/variant.hpp>

#include "types.hpp"


namespace ecell4
{

template <typename T>
struct Quantity
{
    typedef T magnitude_type;
    typedef std::string units_type;

    magnitude_type magnitude;
    units_type units;

    Quantity()
        : magnitude(), units("")
    {
        ;
    }

    Quantity(const magnitude_type& m, const units_type& u="")
        : magnitude(m), units(u)
    {
        ;
    }

    bool operator==(const Quantity& another) const
    {
        return (magnitude == another.magnitude && units == another.units);
    }

    bool operator!=(const Quantity& another) const
    {
        return (magnitude != another.magnitude || units != another.units);
    }
};

template<typename Tstrm_, typename Ttraits_, typename T>
inline std::basic_ostream<Tstrm_, Ttraits_>& operator<<(
    std::basic_ostream<Tstrm_, Ttraits_>& strm, const Quantity<T>& value)
{
    strm << value.magnitude << " " << value.units;
    return strm;
}

} // ecell4

#endif /* ECELL4_QUANTITY_HPP */
