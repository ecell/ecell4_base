#include <algorithm>

#include "exceptions.hpp"
#include "NetworkModel.hpp"


namespace ecell4
{

std::vector<ReactionRule> NetworkModel::query_reaction_rules(
    const Species& sp) const
{
    ReactionRule::reactant_container_type reactants;
    reactants.insert(sp);
    reaction_rules_map_type::const_iterator
        i(reaction_rules_map_.find(reactants));
    std::vector<ReactionRule> retval;
    if (i != reaction_rules_map_.end())
    {
        retval.reserve((*i).second.size());
        for (reaction_rules_map_type::mapped_type::const_iterator
                 j((*i).second.begin()); j != (*i).second.end(); ++j)
        {
            retval.push_back(reaction_rules_[*j]);
        }
    }
    return retval;
}

std::vector<ReactionRule> NetworkModel::query_reaction_rules(
    const Species& sp1, const Species& sp2) const
{
    ReactionRule::reactant_container_type reactants;
    reactants.insert(sp1);
    reactants.insert(sp2);
    reaction_rules_map_type::const_iterator
        i(reaction_rules_map_.find(reactants));
    std::vector<ReactionRule> retval;
    if (i != reaction_rules_map_.end())
    {
        retval.reserve((*i).second.size());
        for (reaction_rules_map_type::mapped_type::const_iterator
                 j((*i).second.begin()); j != (*i).second.end(); ++j)
        {
            retval.push_back(reaction_rules_[*j]);
        }
    }
    return retval;
}

void NetworkModel::add_species(const Species& sp)
{
    if (has_species(sp))
    {
        throw AlreadyExists("species already exists");
    }
    species_.push_back(sp);
}

void NetworkModel::remove_species(const Species& sp)
{
    species_container_type::iterator i(
        std::find(species_.begin(), species_.end(), sp));
    if (i == species_.end())
    {
        throw NotFound("species not found");
    }
    species_.erase(i);
}

bool NetworkModel::has_species(const Species& sp) const
{
    species_container_type::const_iterator i(
        std::find(species_.begin(), species_.end(), sp));
    return (i != species_.end());
}

void NetworkModel::add_reaction_rule(const ReactionRule& rr)
{
    reaction_rule_container_type::const_iterator
        i(std::find(reaction_rules_.begin(), reaction_rules_.end(), rr));
    if (i != reaction_rules_.end())
    {
        throw AlreadyExists("reaction rule already exists");
    }

    reaction_rules_map_[rr.reactants()].insert(reaction_rules_.size());
    reaction_rules_.push_back(rr);
}

void NetworkModel::remove_reaction_rule(const ReactionRule& rr)
{
    reaction_rule_container_type::iterator
        i(std::find(reaction_rules_.begin(), reaction_rules_.end(), rr));
    if (i == reaction_rules_.end())
    {
        throw NotFound("reaction rule not found");
    }

    reaction_rule_container_type::size_type const
        idx(i - reaction_rules_.begin()), last_idx(reaction_rules_.size() - 1);
    reaction_rules_map_type::iterator
        j(reaction_rules_map_.find(rr.reactants()));
    if (j == reaction_rules_map_.end())
    {
        throw IllegalState("no corresponding map key found");
    }
    else if ((*j).second.erase(idx) == 0)
    {
        throw IllegalState("no corresponding map value found");
    }

    if (idx < last_idx)
    {
        reaction_rule_container_type::value_type const
            last_value(reaction_rules_[last_idx]);
        (*i) = last_value;
        j = reaction_rules_map_.find(last_value.reactants());
        if (j == reaction_rules_map_.end())
        {
            throw IllegalState("no corresponding map key for the last found");
        }
        else if ((*j).second.erase(last_idx) == 0)
        {
            throw IllegalState("no corresponding map value for the last found");
        }
        (*j).second.insert(idx);
    }

    reaction_rules_.pop_back();
}

bool NetworkModel::has_reaction_rule(const ReactionRule& rr) const
{
    reaction_rule_container_type::const_iterator
        i(std::find(reaction_rules_.begin(), reaction_rules_.end(), rr));
    return (i != reaction_rules_.end());
}

} // ecell4
