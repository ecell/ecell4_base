#ifndef IDENTIFIER_HPP
#define IDENTIFIER_HPP

#include <utility>

struct default_lot
{
    default_lot& operator=(default_lot const&) {
        return *this;
    }

    bool operator!() const { return true; }

    default_lot& operator++()
    {
        return *this;
    }

    default_lot operator++(int)
    {
        return default_lot();
    }

    default_lot& operator--()
    {
        return *this;
    }

    default_lot operator--(int)
    {
        return default_lot();
    }

    bool operator==(default_lot const& rhs) const
    {
        return true;
    }

    bool operator!=(default_lot const& rhs) const
    {
        return false;
    }

    bool operator<(default_lot const& rhs) const
    {
        return false;
    }

    bool operator>=(default_lot const& rhs) const
    {
        return false;
    }

    bool operator>(default_lot const& rhs) const
    {
        return false;
    }

    bool operator<=(default_lot const& rhs) const
    {
        return false;
    }
};

template<typename Tbase_, typename Tserial_, typename Tlot_ = default_lot>
struct identifier
{
    typedef Tlot_ lot_type;
    typedef Tserial_ serial_type;
    typedef std::pair<lot_type, serial_type> value_type;

    identifier(value_type const& value)
        : value_(value) {}

    Tbase_ lot_add(lot_type const& rhs) const
    {
        return value_type(value_.first + rhs, value_.second);
    }

    Tbase_ lot_subtract(lot_type const& rhs) const
    {
        return value_type(value_.first - rhs, value_.second);
    }

    Tbase_& lot_advance(lot_type const& rhs)
    {
        value_.first += rhs;
        return static_cast<Tbase_&>(*this);
    }

    Tbase_& lot_retrace(lot_type const& rhs)
    {
        value_.first -= rhs;
        return static_cast<Tbase_&>(*this);
    }

    Tbase_ serial_add(serial_type const& rhs) const
    {
        return value_type(value_.first, value_.second + rhs);
    }

    Tbase_ seral_subtract(serial_type const& rhs) const
    {
        return value_type(value_.first, value_.second - rhs);
    }

    Tbase_& serial_advance(serial_type const& rhs)
    {
        value_.second += rhs;
        return static_cast<Tbase_&>(*this);
    }

    Tbase_& serial_retrace(serial_type const& rhs)
    {
        value_.second -= rhs;
        return static_cast<Tbase_&>(*this);
    }

    Tbase_& operator=(Tbase_ const& rhs)
    {
        value_.first = rhs.value_.first;
        value_.second = rhs.value_.second;
    }

    bool operator!() const
    {
        return value_.second == 0;
    }

    bool operator==(Tbase_ const& rhs) const
    {
        return value_.first == rhs.value_.first &&
                value_.second == rhs.value_.second;
    }

    bool operator!=(Tbase_ const& rhs) const
    {
        return value_.first != rhs.value_.first
                || value_.second != rhs.value_.second;
    }

    bool operator<(Tbase_ const& rhs) const
    {
        return value_.second < rhs.value_.second
            || (value_.second == rhs.value_.second &&
                value_.first < rhs.value_.first);
    }

    bool operator>=(Tbase_ const& rhs) const
    {
        return value_.second > rhs.value_.second
            || (value_.second == rhs.value_.second &&
                value_.first >= rhs.value_.first);
    }

    bool operator>(Tbase_ const& rhs) const
    {
        return value_.second > rhs.value_.second
            || (value_.second == rhs.value_.second &&
                value_.first > rhs.value_.first);
    }

    bool operator<=(Tbase_ const& rhs) const
    {
        return value_.second < rhs.value_.second
            || (value_.second == rhs.value_.second &&
                value_.first <= rhs.value_.first);
    }

    operator value_type()
    {
        return value_;
    }

    value_type const& operator()() const
    {
        return value_;
    }

    lot_type& lot()
    {
        return value_.first;
    }

    lot_type const& lot() const
    {
        return value_.first;
    }

    serial_type& serial()
    {
        return value_.second;
    }

    serial_type const& serial() const
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

#endif /* IDENTIFIER_HPP */
