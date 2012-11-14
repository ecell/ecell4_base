#ifndef __MODEL_HPP
#define __MODEL_HPP

#include "types.hpp"

#include "Species.hpp"
#include "ReactionRule.hpp"


namespace ecell4
{

class Model
{
public:

    typedef std::vector<ReactionRule> ReactionRuleVector;

    virtual ReactionRuleVector query_reaction_rules(
        Species const& sp) const = 0;
    virtual ReactionRuleVector query_reaction_rules(
        Species const& sp1, Species const& sp2) const = 0;
};

}

#endif /* __MODEL_HPP */
