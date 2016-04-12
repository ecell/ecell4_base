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
    }
    while (rrexp.next());
    return retval;
}

ReactionRule format_reaction_rule(const ReactionRule& rr)
{
    ReactionRule::reactant_container_type reactants;
    reactants.reserve(rr.reactants().size());
    for (ReactionRule::reactant_container_type::const_iterator i(rr.reactants().begin());
        i != rr.reactants().end(); ++i)
    {
        reactants.push_back(format_species(*i));
    }

    ReactionRule::product_container_type products;
    products.reserve(rr.products().size());
    for (ReactionRule::product_container_type::const_iterator i(rr.products().begin());
        i != rr.products().end(); ++i)
    {
        products.push_back(format_species(*i));
    }

    std::sort(reactants.begin(), reactants.end());
    std::sort(products.begin(), products.end());
    return ReactionRule(reactants, products, rr.k());
    // ReactionRule::reactant_container_type reactants(rr.reactants());
    // ReactionRule::product_container_type products(rr.products());
    // std::sort(reactants.begin(), reactants.end());
    // std::sort(products.begin(), products.end());
    // return ReactionRule(reactants, products, rr.k());
    // return rr;
}

}// ecell4
