#include "Model.hpp"

namespace ecell4
{

ReactionRule create_unimolecular_reaction_rule(
    Species const& reactant1, Species const& product1, double const& k)
{
    ReactionRule rr;
    rr.set_k(k);
    rr.add_reactant(reactant1);
    rr.add_product(product1);
    return rr;
}

ReactionRule create_binding_reaction_rule(
    Species const& reactant1, Species const& reactant2, Species const& product1,
    double const& k)
{
    ReactionRule rr;
    rr.set_k(k);
    rr.add_reactant(reactant1);
    rr.add_reactant(reactant2);
    rr.add_product(product1);
    return rr;
}

ReactionRule create_unbinding_reaction_rule(
    Species const& reactant1, Species const& product1, Species const& product2,
    double const& k)
{
    ReactionRule rr;
    rr.set_k(k);
    rr.add_reactant(reactant1);
    rr.add_product(product1);
    rr.add_product(product2);
    return rr;
}

ReactionRule create_repulsive_reaction_rule(
    Species const& reactant1, Species const& reactant2)
{
    ReactionRule rr;
    rr.set_k(0.0);
    rr.add_reactant(reactant1);
    rr.add_reactant(reactant2);
    return rr;
}

} // ecell4
