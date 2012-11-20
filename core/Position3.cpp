#include "Position3.hpp"


namespace ecell4
{

Position3& Position3::operator+=(Position3 const& rhs)
{
    *this = add(*this, rhs);
    return *this;
}

Position3& Position3::operator-=(Position3 const& rhs)
{
    *this = subtract(*this, rhs);
    return *this;
}

Position3& Position3::operator*=(Position3::value_type const& rhs)
{
    *this = multiply(*this, rhs);
    return *this;
}

Position3& Position3::operator/=(Position3::value_type const& rhs)
{
    *this = divide(*this, rhs);
    return *this;
}

} // ecell4
