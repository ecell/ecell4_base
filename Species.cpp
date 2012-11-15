#include "Species.hpp"


namespace ecell4
{

bool Species::operator==(Species const& rhs) const
{
    return true;
}

bool Species::operator<(Species const& rhs) const
{
    return false;
}

bool Species::operator>(Species const& rhs) const
{
    return false;
}

}
