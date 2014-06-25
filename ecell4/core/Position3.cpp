#include "Position3.hpp"


namespace ecell4
{

Position3& Position3::operator+=(const Position3& rhs)
{
    *this = add(*this, rhs);
    return *this;
}

Position3& Position3::operator-=(const Position3& rhs)
{
    *this = subtract(*this, rhs);
    return *this;
}

Position3& Position3::operator*=(const Position3::value_type& rhs)
{
    *this = multiply(*this, rhs);
    return *this;
}

Position3& Position3::operator/=(const Position3::value_type& rhs)
{
    *this = divide(*this, rhs);
    return *this;
}

} // ecell4
