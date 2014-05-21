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
}


std::string serialize_unit_species_masked(const UnitSpecies& usp)
{
    if (usp.num_sites() == 0)
    {
        return usp.name();
    }

    std::vector<std::string> unstated, stated;
    for (UnitSpecies::container_type::const_iterator i(usp.begin());
        i != usp.end(); ++i)
    {
        const std::string&
            state((*i).second.first), bond((*i).second.second);
        if (state.size() > 0)
        {
            stated.push_back((*i).first + "="
                + (bond.size() > 0? state + "^_" : state));
        }
        else
        {
            unstated.push_back(
                bond.size() > 0? (*i).first + "^_" : (*i).first);
        }
    }

    std::sort(unstated.begin(), unstated.end());
    std::sort(stated.begin(), stated.end());
    return usp.name() + "(" + boost::algorithm::join(unstated, ",")
        + (unstated.size() > 0 && stated.size() > 0? "," : "")
        + boost::algorithm::join(stated, ",") + ")";
}

bool myfunction(const UnitSpecies& lhs, const UnitSpecies& rhs)
{
    // return (lhs < rhs);
    return (serialize_unit_species_masked(lhs) < serialize_unit_species_masked(rhs));
}

std::string serialize_species(const Species& sp)
{
    std::vector<UnitSpecies> units(sp.list_units());
    std::sort(units.begin(), units.end(), myfunction);

    std::vector<UnitSpecies>::const_iterator it(units.begin());
    std::string retval((*it).serial());
    ++it;
    for (; it != units.end(); ++it)
    {
        retval += ".";
        retval += (*it).serial();
    }
    return retval;
}

} // ecell4
