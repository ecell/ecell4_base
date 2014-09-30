#include "Model.hpp"

namespace ecell4
{

ReactionRule create_degradation_reaction_rule(
    const Species& reactant1, const Real& k)
{
    ReactionRule rr;
    rr.set_k(k);
    rr.add_reactant(reactant1);
    return rr;
}

ReactionRule create_synthesis_reaction_rule(
    const Species& product1, const Real& k)
{
    ReactionRule rr;
    rr.set_k(k);
    rr.add_product(product1);
    return rr;
}

ReactionRule create_unimolecular_reaction_rule(
    const Species& reactant1, const Species& product1, const Real& k)
{
    ReactionRule rr;
    rr.set_k(k);
    rr.add_reactant(reactant1);
    rr.add_product(product1);
    return rr;
}

ReactionRule create_binding_reaction_rule(
    const Species& reactant1, const Species& reactant2, const Species& product1,
    const Real& k)
{
    ReactionRule rr;
    rr.set_k(k);
    rr.add_reactant(reactant1);
    rr.add_reactant(reactant2);
    rr.add_product(product1);
    return rr;
}

ReactionRule create_unbinding_reaction_rule(
    const Species& reactant1, const Species& product1, const Species& product2,
    const Real& k)
{
    ReactionRule rr;
    rr.set_k(k);
    rr.add_reactant(reactant1);
    rr.add_product(product1);
    rr.add_product(product2);
    return rr;
}

// ReactionRule create_repulsive_reaction_rule(
//     const Species& reactant1, const Species& reactant2)
// {
//     ReactionRule rr;
//     rr.set_k(0.0);
//     rr.add_reactant(reactant1);
//     rr.add_reactant(reactant2);
//     return rr;
// }

} // ecell4
