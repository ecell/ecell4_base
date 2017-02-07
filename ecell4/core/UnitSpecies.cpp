#include <stdexcept>
#include <boost/algorithm/string.hpp>

#include "UnitSpecies.hpp"

#if defined(HAVE_BOOST_REGEX)
#include <boost/regex.hpp>
#elif defined(WIN32_MSC)
#include <regex>
#else
#include <regex.h>
#endif /* HAVE_BOOST_REGEX */


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
    if (serial == "")
    {
        return;
    }

#if defined(HAVE_BOOST_REGEX) || defined(WIN32_MSC)
#if defined(HAVE_BOOST_REGEX)
    using namespace boost;
#else /* WIN32_MSC */
    using namespace std::tr1;
#endif /* HAVE_BOOST_REGEX */
    regex r1(
        "^\\s*(\\w+)\\s*(\\(\\s*([\\w\\s\\^=,]*)\\))?\\s*$");
    smatch results1;
    if (regex_match(serial, results1, r1))
    {
        name_ = std::string(results1.str(1).c_str());
        if (results1.str(3).size() > 0)
        {
            regex r2(
                "\\s*(\\w+)(\\s*=\\s*(\\w+))?(\\s*\\^\\s*(\\w+))?\\s*");
            // match_results<std::string::const_iterator> results2;
            smatch results2;
            std::vector<std::string> sites;
            boost::split(
                sites, static_cast<const std::string>(results1.str(3)),
                boost::is_any_of(","));
            bool order(false);
            for (std::vector<std::string>::const_iterator i(sites.begin());
                i != sites.end(); ++i)
            {
                if (regex_match(*i, results2, r2))
                {
                    if (results2.str(3).size() > 0)
                    {
                        order = true;
                    }
                    else if (order)
                    {
                        throw std::invalid_argument(
                            "non-keyword arg after keyword arg [" +
                            (*i) + "]"); //XXX:
                    }

                    add_site(
                        results2.str(1), results2.str(3), results2.str(5));
                }
                else
                {
                    throw std::invalid_argument(
                        "a wrong site specification was given [" +
                        (*i) + "]"); //XXX:
                }
            }
        }
    }
    else
    {
        throw std::invalid_argument(
            "a wrong serial was given to UnitSpecies [" + serial + "]"); //XXX:
    }
#else /* regex.h */
    regex_t reg1;
    int errcode = regcomp(&reg1,
        "^[[:blank:]]*([[:alnum:]_]+)[[:blank:]]*"
        "(\\([[:blank:]]*([^\\(\\)[:blank:]][^\\(\\)]*)?\\))?[[:blank:]]*$",
        REG_EXTENDED);
    if (errcode != 0)
    {
        char errbuf[100];
        regerror(errcode, &reg1, errbuf, sizeof(errbuf));
        std::cout << errbuf << std::endl; //XXX: never get here
    }

    regmatch_t match1[4];
    errcode = regexec(&reg1, serial.c_str(), 4, match1, 0);
    if (errcode != 0)
    {
        char errbuf[100];
        regerror(errcode, &reg1, errbuf, sizeof(errbuf));
        throw std::invalid_argument(
            "a wrong serial was given to UnitSpecies [" + serial + "]: "
            + std::string(errbuf)); //XXX:
    }

    name_ = serial.substr(match1[1].rm_so, match1[1].rm_eo - match1[1].rm_so);

    if (match1[3].rm_eo - match1[3].rm_so > 0)
    {
        std::string tmp(
            serial.substr(match1[3].rm_so, match1[3].rm_eo - match1[3].rm_so));

        regex_t reg2;
        errcode = regcomp(&reg2,
            "[[:blank:]]*([[:alnum:]_]+)[[:blank:]]*"
            "(=[[:blank:]]*([[:alnum:]_]+))?[[:blank:]]*"
            "(\\^[[:blank:]]*([[:alnum:]_]+))?[[:blank:]]*(,|$)",
            REG_EXTENDED);
        if (errcode != 0)
        {
            char errbuf[100];
            regerror(errcode, &reg2, errbuf, sizeof(errbuf));
            std::cout << errbuf << std::endl; //XXX: never get here
        }

        regmatch_t match2[7];
        bool order(false);

        while (true)
        {
            errcode = regexec(&reg2, tmp.c_str(), 7, match2, 0);
            if (errcode != 0)
            {
                char errbuf[100];
                regerror(errcode, &reg2, errbuf, sizeof(errbuf));
                throw std::invalid_argument(
                    "wrong site specifiers are given to UnitSpecies ["
                    + serial + "]: " + std::string(errbuf)); //XXX:
            }

            if (match2[3].rm_so != -1)
            {
                order = true;
            }
            else if (order)
            {
                throw std::invalid_argument(
                    "non-keyword arg after keyword arg [" +
                    serial + "]"); //XXX:
            }

            const std::string site_name(
                tmp.substr(match2[1].rm_so, match2[1].rm_eo - match2[1].rm_so));
            const std::string state((match2[3].rm_so != -1)?
                tmp.substr(match2[3].rm_so, match2[3].rm_eo - match2[3].rm_so)
                : "");
            const std::string bond((match2[5].rm_so != -1)?
                tmp.substr(match2[5].rm_so, match2[5].rm_eo - match2[5].rm_so)
                : "");

            add_site(site_name, state, bond);

            if (static_cast<size_t>(match2[0].rm_eo) == tmp.length())
            {
                break;
            }
            else
            {
                tmp = tmp.substr(match2[0].rm_eo);
            }
        }

        regfree(&reg2);
    }

    regfree(&reg1);
#endif /* HAVE_BOOST_REGEX */
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
    return name_ + "(" + boost::algorithm::join(unstated, ",")
        + (unstated.size() > 0 && stated.size() > 0? "," : "")
        + boost::algorithm::join(stated, ",") + ")";

    // std::stringstream unstated, stated;
    // bool is_unstated_empty(true), is_stated_empty(true);
    // unstated << name_ << "(";
    // for (container_type::const_iterator i(sites_.begin());
    //     i != sites_.end(); ++i)
    // {
    //     const std::string& state((*i).second.first);
    //     const std::string& bond((*i).second.second);

    //     if (state.size() > 0)
    //     {
    //         if (is_stated_empty)
    //         {
    //             is_stated_empty = false;
    //         }
    //         else
    //         {
    //             stated << ",";
    //         }
    //         stated << (*i).first << "=" << state;
    //         if (bond.size() > 0)
    //         {
    //             stated << "^" << bond;
    //         }
    //     }
    //     else
    //     {
    //         if (is_unstated_empty)
    //         {
    //             is_unstated_empty = false;
    //         }
    //         else
    //         {
    //             unstated << ",";
    //         }
    //         unstated << (*i).first;
    //         if (bond.size() > 0)
    //         {
    //             unstated << "^" << bond;
    //         }
    //     }
    // }
    // if (!is_unstated_empty && !is_stated_empty)
    // {
    //     unstated << ",";
    // }
    // unstated << stated.str() << ")";
    // return unstated.str();
}

} // ecell4
