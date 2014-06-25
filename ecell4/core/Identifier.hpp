#ifndef __ECELL4_IDENTIFIER_HPP
#define __ECELL4_IDENTIFIER_HPP

#include "config.h"

#include <ostream>
#if defined(HAVE_TR1_FUNCTIONAL)
#include <tr1/functional>
#elif defined(HAVE_STD_HASH)
#include <functional>
#elif defined(HAVE_BOOST_FUNCTIONAL_HASH_HPP)
#include <boost/functional/hash.hpp>
#endif

#include <utility>


namespace ecell4
{

struct DefaultLot
{
    DefaultLot& operator=(const DefaultLot&)
    {
        return *this;
    }

    operator bool() const
    {
        return false;
    }

    bool operator!() const
    {
        return true;
    }

    DefaultLot& operator++()
    {
        return *this;
    }

    DefaultLot operator++(int)
    {
        return DefaultLot();
    }

    DefaultLot& operator--()
    {
        return *this;
    }

    DefaultLot operator--(int)
    {
        return DefaultLot();
    }

    bool operator==(const DefaultLot& rhs) const
    {
        return true;
    }

    bool operator!=(const DefaultLot& rhs) const
    {
        return false;
    }

    bool operator<(const DefaultLot& rhs) const
    {
        return false;
    }

    bool operator>=(const DefaultLot& rhs) const
    {
        return false;
    }

    bool operator>(const DefaultLot& rhs) const
    {
        return false;
    }

    bool operator<=(const DefaultLot& rhs) const
    {
        return false;
    }
};

template<typename Tbase_, typename Tserial_, typename Tlot_ = DefaultLot>
struct Identifier
{
public:

    typedef Tlot_ lot_type;
    typedef Tserial_ serial_type;
    typedef std::pair<lot_type, serial_type> value_type;

public:

    Identifier(const value_type& value)
        : value_(value)
    {
        ;
    }

    Tbase_ lot_add(const lot_type& rhs) const
    {
        return value_type(value_.first + rhs, value_.second);
    }

    Tbase_ lot_subtract(const lot_type& rhs) const
    {
        return value_type(value_.first - rhs, value_.second);
    }

    Tbase_& lot_advance(const lot_type& rhs)
    {
        value_.first += rhs;
        return static_cast<Tbase_&>(*this);
    }

    Tbase_& lot_retrace(const lot_type& rhs)
    {
        value_.first -= rhs;
        return static_cast<Tbase_&>(*this);
    }

    Tbase_ serial_add(const serial_type& rhs) const
    {
        return value_type(value_.first, value_.second + rhs);
    }

    Tbase_ seral_subtract(const serial_type& rhs) const
    {
        return value_type(value_.first, value_.second - rhs);
    }

    Tbase_& serial_advance(const serial_type& rhs)
    {
        value_.second += rhs;
        return static_cast<Tbase_&>(*this);
    }

    Tbase_& serial_retrace(const serial_type& rhs)
    {
        value_.second -= rhs;
        return static_cast<Tbase_&>(*this);
    }

    Tbase_& operator=(const Tbase_& rhs)
    {
        value_.first = rhs.value_.first;
        value_.second = rhs.value_.second;
    }

    operator bool() const
    {
        return value_.second != 0;
    }

    bool operator!() const
    {
        return value_.second == 0;
    }

    bool operator==(const Tbase_& rhs) const
    {
        return value_.first == rhs.value_.first &&
            value_.second == rhs.value_.second;
    }

    bool operator!=(const Tbase_& rhs) const
    {
        return value_.first != rhs.value_.first
            || value_.second != rhs.value_.second;
    }

    bool operator<(const Tbase_& rhs) const
    {
        return value_.second < rhs.value_.second
            || (value_.second == rhs.value_.second &&
                value_.first < rhs.value_.first);
    }

    bool operator>=(const Tbase_& rhs) const
    {
        return value_.second > rhs.value_.second
            || (value_.second == rhs.value_.second &&
                value_.first >= rhs.value_.first);
    }

    bool operator>(const Tbase_& rhs) const
    {
        return value_.second > rhs.value_.second
            || (value_.second == rhs.value_.second &&
                value_.first > rhs.value_.first);
    }

    bool operator<=(const Tbase_& rhs) const
    {
        return value_.second < rhs.value_.second
            || (value_.second == rhs.value_.second &&
                value_.first <= rhs.value_.first);
    }

    operator value_type() const
    {
        return value_;
    }

    const value_type& operator()() const
    {
        return value_;
    }

    lot_type& lot()
    {
        return value_.first;
    }

    const lot_type& lot() const
    {
        return value_.first;
    }

    serial_type& serial()
    {
        return value_.second;
    }

    const serial_type& serial() const
    {
        return value_.second;
    }

    value_type& operator()()
    {
        return value_;
    }

protected:

    value_type value_;
};

struct ParticleID:
        public Identifier<ParticleID, unsigned long long, int>
{
    typedef Identifier<ParticleID, unsigned long long, int> base_type;

    ParticleID(const value_type& value = value_type(0, 0))
        : base_type(value)
    {
        ;
    }
};

} // ecell4

#if defined(HAVE_TR1_FUNCTIONAL)
namespace std
{

namespace tr1
{
#elif defined(HAVE_STD_HASH)
namespace std
{
#elif defined(HAVE_BOOST_FUNCTIONAL_HASH_HPP)
namespace boost
{
#endif

template<>
struct hash<ecell4::ParticleID>
{
    std::size_t operator()(const ecell4::ParticleID& val) const
    {
        return static_cast<std::size_t>(val().first ^ val().second);
    }
};

#if defined(HAVE_TR1_FUNCTIONAL)
} // tr1

} // std
#elif defined(HAVE_STD_HASH)
} // std
#elif defined(HAVE_BOOST_FUNCTIONAL_HASH_HPP)
} // boost
#endif

#endif /* __ECELL4_IDENTIFIER_HPP */
