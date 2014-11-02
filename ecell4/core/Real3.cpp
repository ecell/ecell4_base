#include "Real3.hpp"


namespace ecell4
{

Real3& Real3::operator+=(const Real3& rhs)
{
    *this = add(*this, rhs);
    return *this;
}

Real3& Real3::operator-=(const Real3& rhs)
{
    *this = subtract(*this, rhs);
    return *this;
}

Real3& Real3::operator*=(const Real3::value_type& rhs)
{
    *this = multiply(*this, rhs);
    return *this;
}

Real3& Real3::operator/=(const Real3::value_type& rhs)
{
    *this = divide(*this, rhs);
    return *this;
}

} // ecell4
