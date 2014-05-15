#include "Context.hpp"
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

bool Species::match(const Species& target) const
{
    return spmatch(*this, target);
    // container_type::const_iterator i(units_.begin()), j(target.begin());
    // while (i != units_.end())
    // {
    //     const UnitSpecies& usp(*i);

    //     j = std::lower_bound(j, target.end(), usp);
    //     if (j == target.end())
    //     {
    //         return false;
    //     }

    //     const container_type::const_iterator
    //         nexti(std::upper_bound(i, units_.end(), usp)),
    //         nextj(std::upper_bound(j, target.end(), usp));
    //     if (nextj - j < nexti - i)
    //     {
    //         return false;
    //     }

    //     i = nexti;
    //     j = nextj;
    // }
    // return true;
}

} // ecell4
