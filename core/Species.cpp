#include "Species.hpp"


namespace ecell4
{

bool Species::operator==(Species const& rhs) const
{
    return (name() == rhs.name());
}

bool Species::operator<(Species const& rhs) const
{
    return (name() < rhs.name());
}

bool Species::operator>(Species const& rhs) const
{
    return (name() > rhs.name());
}

} // ecell4
