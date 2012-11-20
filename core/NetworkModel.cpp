#include "NetworkModel.hpp"


namespace ecell4
{

ReactionRuleVector NetworkModel::query_reaction_rules(
    Species const& sp) const
{
    ReactionRuleVector retval;
    return retval;
}

ReactionRuleVector NetworkModel::query_reaction_rules(
    Species const& sp1, Species const& sp2) const
{
    ReactionRuleVector retval;
    return retval;
}

bool NetworkModel::add_species(Species const& sp)
{
    return true;
}

bool NetworkModel::add_reaction_rule(ReactionRule const& rr)
{
    return true;
}

} // ecell4
