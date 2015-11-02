#include <sstream>
#include <boost/algorithm/string.hpp>

#include "ReactionRule.hpp"
#include "Context.hpp"

namespace ecell4
{

const std::string ReactionRule::as_string() const
{
    std::stringstream oss;
    std::vector<std::string> tmp;
    for (reactant_container_type::const_iterator i(reactants_.begin());
        i != reactants_.end(); ++i)
    {
        tmp.push_back((*i).serial());
    }
    oss << boost::algorithm::join(tmp, "+") << ">";
    tmp.clear();
    for (product_container_type::const_iterator i(products_.begin());
        i != products_.end(); ++i)
    {
        tmp.push_back((*i).serial());
    }
    oss << boost::algorithm::join(tmp, "+") << "|" << k_;
    return oss.str();
}

Integer ReactionRule::count(const reactant_container_type& reactants) const
{
    ReactionRuleExpressionMatcher rrexp(*this);
    if (!rrexp.match(reactants))
    {
        return 0;
    }
    Integer n(1);
    while (rrexp.next())
    {
        ++n;
    }
    return n;
}

std::vector<ReactionRule> ReactionRule::generate(const reactant_container_type& reactants) const
{
    ReactionRuleExpressionMatcher rrexp(*this);
    std::vector<ReactionRule> retval;
    if (!rrexp.match(reactants))
    {
        return retval;
    }

    do
    {
        const ReactionRule rr(reactants, rrexp.generate(), this->k());
        std::vector<ReactionRule>::iterator
            i(std::find(retval.begin(), retval.end(), rr));
        if (i != retval.end())
        {
            ;
        }
        else
        {
            retval.push_back(rr);
        }
        // if (i != retval.end())
        // {
        //     (*i).set_k((*i).k() + rr.k());
        // }
        // else
        // {
        //     retval.push_back(rr);
        // }
    }
    while (rrexp.next());
    return retval;

    // const std::vector<std::vector<Species> > possibles(rrgenerate(*this, reactants));
    // std::vector<ReactionRule> retval;
    // retval.reserve(possibles.size());
    // for (std::vector<std::vector<Species> >::const_iterator i(possibles.begin());
    //     i != possibles.end(); ++i)
    // {
    //     const ReactionRule rr(reactants, *i, this->k());
    //     std::vector<ReactionRule>::iterator
    //         j(std::find(retval.begin(), retval.end(), rr));
    //     if (j != retval.end())
    //     {
    //         (*j).set_k((*j).k() + rr.k());
    //     }
    //     else
    //     {
    //         retval.push_back(rr);
    //     }
    // }
    // return retval;
}

}// ecell4
