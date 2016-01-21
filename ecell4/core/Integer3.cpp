#include "Integer3.hpp"

namespace ecell4 {

Integer3& Integer3::operator+=(const Integer3& rhs)
{
    *this = add(*this, rhs);
    return *this;
}

Integer3& Integer3::operator-=(const Integer3& rhs)
{
    *this = subtract(*this, rhs);
    return *this;
}

Integer3& Integer3::operator*=(const Integer3::value_type& rhs)
{
    *this = multiply(*this, rhs);
    return *this;
}

Integer3 Integer3::east() const
{
    Integer3 retval(*this);
    retval.col += 1;
    return retval;
}

Integer3 Integer3::west() const
{
    Integer3 retval(*this);
    retval.col -= 1;
    return retval;
}

Integer3 Integer3::south() const
{
    Integer3 retval(*this);
    retval.row += 1;
    return retval;
}

Integer3 Integer3::north() const
{
    Integer3 retval(*this);
    retval.row -= 1;
    return retval;
}

Integer3 Integer3::dorsal() const
{
    Integer3 retval(*this);
    retval.layer += 1;
    return retval;
}

Integer3 Integer3::ventral() const
{
    Integer3 retval(*this);
    retval.layer -= 1;
    return retval;
}

} // ecell4
