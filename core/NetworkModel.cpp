#include <algorithm>

#include "exceptions.hpp"
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
    if (has_species(sp))
    {
        throw AlreadyExists("species already exists");
    }
    species_.push_back(sp);
    return true;
}

bool NetworkModel::has_species(Species const& sp) const
{
    typename species_container_type::const_iterator i(
        std::find(species_.begin(), species_.end(), sp));
    return (i != species_.end());
}

bool NetworkModel::add_reaction_rule(ReactionRule const& rr)
{
    std::pair<typename reaction_rules_type::mapped_type::iterator, bool>
        retval(reaction_rules_[rr.reactants()].insert(rr));
    if (!retval.second)
    {
        throw AlreadyExists("reaction rule already exists");
    }
    return true;
}

void NetworkModel::remove_reaction_rule(ReactionRule const& rr)
{
    typename reaction_rules_type::iterator
        i(reaction_rules_.find(rr.reactants()));
    if (i == reaction_rules_.end())
    {
        throw NotFound("reaction rule not found");
    }
    else if ((*i).second.erase(rr) == 0) /// remove a reaction rule here
    {
        throw NotFound("reaction rule not found");
    }
}

} // ecell4
