#ifndef __MODEL_HPP
#define __MODEL_HPP

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

    typedef Model::ReactionRuleVector ReactionRuleVector;

    virtual bool add_species(Species const& sp) = 0;
    virtual bool add_reaction_rule(ReactionRule const& rr) = 0;
};

}

#endif /* __MODEL_HPP */
