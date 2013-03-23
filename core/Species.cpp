#include "Species.hpp"


namespace ecell4
{

bool Species::operator==(const Species& rhs) const
{
    return (serial() == rhs.serial());
}

bool Species::operator<(const Species& rhs) const
{
    return (serial() < rhs.serial());
}

bool Species::operator>(const Species& rhs) const
{
    return (serial() > rhs.serial());
}

} // ecell4
