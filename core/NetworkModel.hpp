#ifndef __NETWORK_MODEL_HPP
#define __NETWORK_MODEL_HPP

#include <set>

#include "types.hpp"

#include "Species.hpp"
#include "ReactionRule.hpp"
#include "Model.hpp"


namespace ecell4
{

class NetworkModel
    : public Model
{
public:

    typedef std::vector<Species> species_container_type;

    typedef
    std::map<ReactionRule::reactants_type, std::set<ReactionRule> >
    reaction_rules_type;

    NetworkModel()
        : species_(), reaction_rules_()
    {
        ;
    }

    ReactionRuleVector query_reaction_rules(Species const& sp) const;
    ReactionRuleVector query_reaction_rules(
        Species const& sp1, Species const& sp2) const;

    bool add_species(Species const& sp);
    void remove_species(Species const& sp);
    bool has_species(Species const& sp) const;

    bool add_reaction_rule(ReactionRule const& rr);
    void remove_reaction_rule(ReactionRule const& rr);
    bool has_reaction_rule(ReactionRule const& rr) const;

protected:

    species_container_type species_;
    reaction_rules_type reaction_rules_;
};

} // ecell4

#endif /* __NETWORK_MODEL_HPP */
