#include <boost/algorithm/string.hpp>
#include <boost/regex.hpp>

#include "UnitSpecies.hpp"


namespace ecell4
{

void UnitSpecies::clear()
{
    name_ = "";
    sites_.clear();
}

void UnitSpecies::deserialize(const UnitSpecies::serial_type& serial)
{
    clear();

    boost::regex r1(
        "^\\s*(\\w+)\\s*(\\(\\s*([\\w\\s\\^=,]*)\\))?\\s*$");
    boost::match_results<std::string::const_iterator> results1;
    if (boost::regex_match(serial, results1, r1))
    {
        name_ = results1.str(1);
        if (results1.str(3).size() > 0)
        {
            boost::regex r2(
                "\\s*(\\w+)(\\s*=\\s*(\\w+))?(\\s*\\^\\s*(\\w+))?\\s*");
            boost::match_results<std::string::const_iterator> results2;
            std::vector<std::string> sites;
            boost::split(
                sites, static_cast<const std::string>(results1.str(3)),
                boost::is_any_of(","));
            bool order(false);
            for (std::vector<std::string>::const_iterator i(sites.begin());
                i != sites.end(); ++i)
            {
                if (boost::regex_match(*i, results2, r2))
                {
                    if (results2.str(3).size() > 0)
                    {
                        order = true;
                    }
                    else if (order)
                    {
                        //XXX:
                    }

                    add_site(
                        results2.str(1), results2.str(3), results2.str(5));
                }
                else
                {
                    //XXX:
                }
            }
        }
    }
    else
    {
        name_ = "FAIL"; //XXX:
    }
}

UnitSpecies::serial_type UnitSpecies::serial() const
{
    if (sites_.size() == 0)
    {
        return name_;
    }

    std::vector<std::string> unstated, stated;
    for (container_type::const_iterator i(sites_.begin());
        i != sites_.end(); ++i)
    {
        const std::string&
            state((*i).second.first), bond((*i).second.second);
        if (state.size() > 0)
        {
            stated.push_back((*i).first + "="
                + (bond.size() > 0? state + "^" + bond : state));
        }
        else
        {
            unstated.push_back(
                bond.size() > 0? (*i).first + "^" + bond : (*i).first);
        }
    }

    std::sort(unstated.begin(), unstated.end());
    std::sort(stated.begin(), stated.end());
    return name_ + "(" + boost::algorithm::join(unstated, ",")
        + (unstated.size() > 0 && stated.size() > 0? "," : "")
        + boost::algorithm::join(stated, ",") + ")";
}

} // ecell4
