#include <algorithm>

#include "exceptions.hpp"
#include "NetworkModel.hpp"


namespace ecell4
{

ReactionRuleVector NetworkModel::query_reaction_rules(Species const& sp) const
{
    ReactionRule::reactants_type reactants;
    reactants.insert(sp);
    reaction_rules_type::const_iterator i(reaction_rules_.find(reactants));
    ReactionRuleVector retval((*i).second.begin(), (*i).second.end());
    return retval;
}

ReactionRuleVector NetworkModel::query_reaction_rules(
    Species const& sp1, Species const& sp2) const
{
    ReactionRule::reactants_type reactants;
    reactants.insert(sp1);
    reactants.insert(sp2);
    reaction_rules_type::const_iterator i(reaction_rules_.find(reactants));
    ReactionRuleVector retval((*i).second.begin(), (*i).second.end());
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

void NetworkModel::remove_species(Species const& sp)
{
    species_container_type::iterator i(
        std::find(species_.begin(), species_.end(), sp));
    if (i == species_.end())
    {
        throw NotFound("species not found");
    }
    species_.erase(i);
}

bool NetworkModel::has_species(Species const& sp) const
{
    species_container_type::const_iterator i(
        std::find(species_.begin(), species_.end(), sp));
    return (i != species_.end());
}

bool NetworkModel::add_reaction_rule(ReactionRule const& rr)
{
    std::pair<reaction_rules_type::mapped_type::iterator, bool>
        retval(reaction_rules_[rr.reactants()].insert(rr));
    if (!retval.second)
    {
        throw AlreadyExists("reaction rule already exists");
    }
    return true;
}

void NetworkModel::remove_reaction_rule(ReactionRule const& rr)
{
    reaction_rules_type::iterator
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

bool NetworkModel::has_reaction_rule(ReactionRule const& rr) const
{
    reaction_rules_type::const_iterator
        i(reaction_rules_.find(rr.reactants()));
    if (i == reaction_rules_.end())
    {
        return false;
    }
    return ((*i).second.find(rr) != (*i).second.end());
}

} // ecell4
