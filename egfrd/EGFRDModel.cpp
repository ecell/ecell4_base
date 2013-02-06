#include <algorithm>

#include "NetworkRules.hpp"
#include "ReactionRule.hpp"

#include <ecell4/core/exceptions.hpp>

#include "EGFRDModel.hpp"


namespace ecell4
{

namespace egfrd
{

std::vector<ReactionRule> EGFRDModel::query_reaction_rules(
    Species const& sp) const
{
    return network_model_.query_reaction_rules(sp);
}

std::vector<ReactionRule> EGFRDModel::query_reaction_rules(
    Species const& sp1, Species const& sp2) const
{
    return network_model_.query_reaction_rules(sp1, sp2);
}

void EGFRDModel::add_species(Species const& sp)
{
    network_model_.add_species(sp);

    boost::shared_ptr< ::SpeciesType> st(new ::SpeciesType());
    (*st)["name"] = boost::lexical_cast<std::string>(sp.name());
    (*st)["D"] = boost::lexical_cast<std::string>(sp.get_attribute("D"));
    (*st)["radius"] = boost::lexical_cast<std::string>(
        sp.get_attribute("radius"));
    particle_model_.add_species_type(st);

    sid_map_[sp.serial()] = st->id();
}

void EGFRDModel::remove_species(Species const& sp)
{
    throw NotImplemented("not implemented yet.");

    // network_model_.remove_species(sp);
}

bool EGFRDModel::has_species(Species const& sp) const
{
    return network_model_.has_species(sp);
}

void EGFRDModel::add_reaction_rule(ReactionRule const& rr)
{
    network_model_.add_reaction_rule(rr);

    std::vector< ::SpeciesTypeID> products;
    for (ReactionRule::reactant_container_type::const_iterator
             j(rr.products().begin()); j != rr.products().end(); ++j)
    {
        products.push_back(get_species_type_id(*j));
    }

    ReactionRule::reactant_container_type::const_iterator
        r(rr.reactants().begin());
    switch (rr.reactants().size())
    {
    case 1:
        {
            ::SpeciesTypeID const sid1(get_species_type_id(*r));
            particle_model_.network_rules().add_reaction_rule(
                ::new_reaction_rule(sid1, products, rr.k()));
        }
        break;
    case 2:
        {
            ::SpeciesTypeID const sid1(get_species_type_id(*r));
            ++r;
            ::SpeciesTypeID const sid2(get_species_type_id(*r));
            particle_model_.network_rules().add_reaction_rule(
                ::new_reaction_rule(sid1, sid2, products, rr.k()));
        }
        break;
    default:
        throw NotSupported("the number of reactants must be 1 or 2.");
    }
}

void EGFRDModel::remove_reaction_rule(ReactionRule const& rr)
{
    throw NotImplemented("not implemented yet.");

    // network_model_.remove_reaction_rule(rr);
}

bool EGFRDModel::has_reaction_rule(ReactionRule const& rr) const
{
    return network_model_.has_reaction_rule(rr);
}

} // egfrd

} // ecell4
