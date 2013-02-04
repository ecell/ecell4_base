#include "Model.hpp"

namespace ecell4
{

ReactionRule create_association_reaction_rule(
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

ReactionRule create_dissociation_reaction_rule(
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

} // ecell4
