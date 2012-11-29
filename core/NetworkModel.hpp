#ifndef __NETWORK_MODEL_HPP
#define __NETWORK_MODEL_HPP

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

    NetworkModel()
    {
        ;
    }

    ReactionRuleVector query_reaction_rules(Species const& sp) const;
    ReactionRuleVector query_reaction_rules(
        Species const& sp1, Species const& sp2) const;

    bool add_species(Species const& sp);
    bool has_species(Species const& sp) const;
    bool add_reaction_rule(ReactionRule const& rr);

protected:

    SpeciesVector species_;
};

} // ecell4

#endif /* __NETWORK_MODEL_HPP */
