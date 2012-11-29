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
    SpeciesVector::const_iterator i(
        std::find(species_.begin(), species_.end(), sp));
    return (i != species_.end());
}

bool NetworkModel::add_reaction_rule(ReactionRule const& rr)
{
    return true;
}

} // ecell4
