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
    return ReactionRuleExpressionMatcher(*this).gen(reactants);
}

ReactionRule format_reaction_rule_with_nosort(const ReactionRule& rr)
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

    return ReactionRule(reactants, products, rr.k());
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

ReactionRule create_degradation_reaction_rule(
    const Species& reactant1, const Real& k)
{
    ReactionRule rr;
    rr.set_k(k);
    rr.add_reactant(reactant1);
    return rr;
}

ReactionRule create_synthesis_reaction_rule(
    const Species& product1, const Real& k)
{
    ReactionRule rr;
    rr.set_k(k);
    rr.add_product(product1);
    return rr;
}

ReactionRule create_unimolecular_reaction_rule(
    const Species& reactant1, const Species& product1, const Real& k)
{
    ReactionRule rr;
    rr.set_k(k);
    rr.add_reactant(reactant1);
    rr.add_product(product1);
    return rr;
}

ReactionRule create_binding_reaction_rule(
    const Species& reactant1, const Species& reactant2, const Species& product1,
    const Real& k)
{
    ReactionRule rr;
    rr.set_k(k);
    rr.add_reactant(reactant1);
    rr.add_reactant(reactant2);
    rr.add_product(product1);
    return rr;
}

ReactionRule create_unbinding_reaction_rule(
    const Species& reactant1, const Species& product1, const Species& product2,
    const Real& k)
{
    ReactionRule rr;
    rr.set_k(k);
    rr.add_reactant(reactant1);
    rr.add_product(product1);
    rr.add_product(product2);
    return rr;
}

// ReactionRule create_repulsive_reaction_rule(
//     const Species& reactant1, const Species& reactant2)
// {
//     ReactionRule rr;
//     rr.set_k(0.0);
//     rr.add_reactant(reactant1);
//     rr.add_reactant(reactant2);
//     return rr;
// }

}// ecell4
