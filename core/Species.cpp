#include "Species.hpp"


namespace ecell4
{

bool Species::operator==(Species const& rhs) const
{
    return (serial() == rhs.serial());
}

bool Species::operator<(Species const& rhs) const
{
    return (serial() < rhs.serial());
}

bool Species::operator>(Species const& rhs) const
{
    return (serial() > rhs.serial());
}

} // ecell4
