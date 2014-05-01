#ifndef __ECELL4_MODEL_HPP
#define __ECELL4_MODEL_HPP

#include "types.hpp"
#include "Species.hpp"
#include "ReactionRule.hpp"
#include "exceptions.hpp"


namespace ecell4
{

ReactionRule create_unimolecular_reaction_rule(
    const Species& reactant1, const Species& product1, const Real& k);

ReactionRule create_binding_reaction_rule(
    const Species& reactant1, const Species& reactant2, const Species& product1,
    const Real& k);

ReactionRule create_unbinding_reaction_rule(
    const Species& reactant1, const Species& product1, const Species& product2,
    const Real& k);

ReactionRule create_degradation_reaction_rule(
    const Species& reactant1, const Real& k);

ReactionRule create_synthesis_reaction_rule(
    const Species& product1, const Real& k);

// ReactionRule create_repulsive_reaction_rule(
//     const Species& reactant1, const Species& reactant2);

class Model
{
public:

    virtual ~Model()
    {
        ;
    }

    // ModelTraits

    /**
     * a fundamental function to query unimolecular reaction rules from a reactant.
     * this must be overloaded by any sub classes of Model.
     * @param species Species of a reactant
     * @return the vector of ReactionRule(s)
     */
    virtual std::vector<ReactionRule> query_reaction_rules(
        const Species& sp) const = 0;

    /**
     * a fundamental function to query bimolecular reaction rules from reactants.
     * this must be overloaded by any sub classes of Model.
     * @param species1 Species of the first reactant
     * @param species2 Species of the second reactant
     * @return the vector of ReactionRule(s)
     */
    virtual std::vector<ReactionRule> query_reaction_rules(
        const Species& sp1, const Species& sp2) const = 0;

    // NetworkModelTraits

    /**
     * add attributes of species to the model.
     * this function is a part of the trait of Model.
     * @param species a new Species
     */
    virtual void add_species_attribute(const Species& sp)
    {
        throw NotSupported(
            "add_species_attribute is not supported in this model class");
    }

    /**
     * return if a species attribute is in the model, or not.
     * this function is a part of the trait of Model.
     * @param species a Species
     * @return if the species exists, or not
     */
    virtual bool has_species_attribute(const Species& sp) const
    {
        throw NotSupported(
            "has_species_attribute is not supported in this model class");
    }

    /**
     * remove attributes of species in the model.
     * this function is a part of the trait of Model.
     * @param species a new Species
     */
    virtual void remove_species_attribute(const Species& sp)
    {
        throw NotSupported(
            "remove_species_attribute is not supported in this model class");
    }

    /**
     * apply attributes of species to the given species.
     * this function is a part of the trait of Model.
     * @param species an original Species
     */
    virtual Species apply_species_attributes(const Species& sp) const
    {
        throw NotSupported(
            "apply_species_attributes is not supported in this model class");
    }

    /**
     * add a reaction rule to the model.
     * this function is a part of the trait of Model.
     * @param rr a new ReactionRule
     * @return if the reaction rule is not registered yet.
     */
    virtual void add_reaction_rule(const ReactionRule& rr)
    {
        throw NotSupported(
            "add_reaction_rule is not supported in this model class");
    }

    /**
     * remove a reaction rule in the model.
     * this function is a part of the trait of Model.
     * @param rr a new ReactionRule
     */
    virtual void remove_reaction_rule(const ReactionRule& rr)
    {
        throw NotSupported(
            "remove_reaction_rule is not supported in this model class");
    }

    /**
     * return if a reaction rule is in the model, or not.
     * this function is a part of the trait of Model.
     * @param rr a reaction rule
     * @return if the reaction rule exists, or not
     */
    virtual bool has_reaction_rule(const ReactionRule& rr)
    {
        throw NotSupported(
            "has_reaction_rule is not supported in this model class");
    }
};

} // ecell4

#endif /* __ECELL4_MODEL_HPP */
