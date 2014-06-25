#include <algorithm>

#include "exceptions.hpp"
#include "NetworkModel.hpp"


namespace ecell4
{

std::vector<ReactionRule> NetworkModel::query_reaction_rules(
    const Species& sp) const
{
    ReactionRule::reactant_container_type reactants;
    reactants.push_back(sp);
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
    std::vector<ReactionRule> retval;
    ReactionRule::reactant_container_type reactants;
    reactants.push_back(sp1);
    reactants.push_back(sp2);

    reaction_rules_map_type::const_iterator
        i(reaction_rules_map_.find(reactants));
    if (i != reaction_rules_map_.end())
    {
        retval.reserve((*i).second.size());
        for (reaction_rules_map_type::mapped_type::const_iterator
                 j((*i).second.begin()); j != (*i).second.end(); ++j)
        {
            retval.push_back(reaction_rules_[*j]);
        }
    }

    std::swap(reactants[0], reactants[1]);
    i = reaction_rules_map_.find(reactants);
    if (i != reaction_rules_map_.end())
    {
        retval.reserve(retval.size() + (*i).second.size());
        for (reaction_rules_map_type::mapped_type::const_iterator
                 j((*i).second.begin()); j != (*i).second.end(); ++j)
        {
            retval.push_back(reaction_rules_[*j]);
        }
    }
    return retval;
}

void NetworkModel::initialize()
{
    if (!dirty_)
    {
        return; // do nothing
    }

    species_cache_.clear();
    for (reaction_rule_container_type::const_iterator
        i(reaction_rules_.begin()); i != reaction_rules_.end(); ++i)
    {
        const ReactionRule::reactant_container_type&
            reactants((*i).reactants());
        const ReactionRule::product_container_type&
            products((*i).products());
        std::copy(reactants.begin(), reactants.end(),
                  std::back_inserter(species_cache_));
        std::copy(products.begin(), products.end(),
                  std::back_inserter(species_cache_));
    }
    std::sort(species_cache_.begin(), species_cache_.end());
    species_cache_.erase(
        std::unique(species_cache_.begin(), species_cache_.end()),
        species_cache_.end());

    dirty_ = false;
}

void NetworkModel::add_species_attribute(const Species& sp)
{
    if (has_species_attribute(sp))
    {
        throw AlreadyExists("species already exists");
    }
    species_attributes_.push_back(sp);

    dirty_ = true;
}

void NetworkModel::remove_species_attribute(const Species& sp)
{
    species_container_type::iterator i(
        std::find(species_attributes_.begin(), species_attributes_.end(), sp));
    if (i == species_attributes_.end())
    {
        std::ostringstream message;
        message << "Speices [" << sp.serial() << "] not found";
        throw NotFound(message.str()); // use boost::format if it's allowed
    }
    species_attributes_.erase(i);

    dirty_ = true;
}

bool NetworkModel::has_species_attribute(const Species& sp) const
{
    species_container_type::const_iterator i(
        std::find(species_attributes_.begin(), species_attributes_.end(), sp));
    return (i != species_attributes_.end());
}

void NetworkModel::add_reaction_rule(const ReactionRule& rr)
{
    reaction_rule_container_type::const_iterator
        i(std::find(reaction_rules_.begin(), reaction_rules_.end(), rr));
    if (i != reaction_rules_.end())
    {
        throw AlreadyExists("reaction rule already exists");
    }

    const reaction_rule_container_type::size_type idx(reaction_rules_.size());
    reaction_rules_map_[rr.reactants()].insert(idx);
    reaction_rules_.push_back(rr);

    dirty_ = true;
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
    dirty_ = true;
}

bool NetworkModel::has_reaction_rule(const ReactionRule& rr) const
{
    reaction_rule_container_type::const_iterator
        i(std::find(reaction_rules_.begin(), reaction_rules_.end(), rr));
    return (i != reaction_rules_.end());
}

} // ecell4
