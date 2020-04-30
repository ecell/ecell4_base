#ifndef ECELL4_NETFREE_MODEL_HPP
#define ECELL4_NETFREE_MODEL_HPP

#include <map>
#include <set>
#include <algorithm>
#include <iterator>
#include <memory>

#include "types.hpp"
#include "Species.hpp"
#include "ReactionRule.hpp"
#include "Model.hpp"

#include "Context.hpp"
#include "NetworkModel.hpp"


namespace ecell4
{

class NetfreeModel
    : public Model
{
public:

    typedef Model base_type;
    typedef base_type::species_container_type species_container_type;
    typedef base_type::reaction_rule_container_type reaction_rule_container_type;

public:

    NetfreeModel()
        : base_type(), species_attributes_(), species_attributes_proceed_(), reaction_rules_(), effective_(false)
    {
        ;
    }

    virtual ~NetfreeModel()
    {
        ;
    }

    // ModelTraits

    std::vector<ReactionRule> query_reaction_rules(const Species& sp) const;
    std::vector<ReactionRule> query_reaction_rules(
        const Species& sp1, const Species& sp2) const;

    std::vector<ReactionRule> query_reaction_rules(
        const std::vector<Species>& splist, const std::vector<Species>::size_type n) const;

    inline std::vector<ReactionRule> query_reaction_rules(
        const std::vector<Species>& reactants) const
    {
        return this->query_reaction_rules(reactants, reactants.size());
    }

    Integer apply(const Species& pttrn, const Species& sp) const;
    std::vector<ReactionRule> apply(
        const ReactionRule& rr,
        const ReactionRule::reactant_container_type& reactants) const;

    Species apply_species_attributes(const Species& sp) const;

    // NetfreeModelTraits

    bool update_species_attribute(const Species& sp);
    void add_species_attribute(const Species& sp, const bool proceed = false);
    bool has_species_attribute(const Species& sp) const;
    bool has_species_attribute_exact(const Species& sp) const;
    void remove_species_attribute(const Species& sp);

    void add_reaction_rule(const ReactionRule& rr);
    void remove_reaction_rule(const ReactionRule& rr);
    bool has_reaction_rule(const ReactionRule& rr) const;

    const species_container_type& species_attributes() const
    {
        return species_attributes_;
    }

    const std::vector<bool>& species_attributes_proceed() const
    {
        return species_attributes_proceed_;
    }

    const reaction_rule_container_type& reaction_rules() const
    {
        return reaction_rules_;
    }

    // Optional functions

    std::shared_ptr<Model> expand(
        const std::vector<Species>& sp, const Integer max_itr,
        const std::map<Species, Integer>& max_stoich) const;
    std::shared_ptr<Model> expand(
        const std::vector<Species>& sp, const Integer max_itr) const;
    std::shared_ptr<Model> expand(const std::vector<Species>& sp) const;

    void set_effective(const bool effective)
    {
        effective_ = effective;
    }

    const bool effective() const
    {
        return effective_;
    }

protected:

    species_container_type species_attributes_;
    std::vector<bool> species_attributes_proceed_;  //XXX:
    reaction_rule_container_type reaction_rules_;

    bool effective_;
};

namespace extras
{

std::pair<std::shared_ptr<NetworkModel>, bool> generate_network_from_netfree_model(
    const NetfreeModel& nfm, const std::vector<Species>& seeds, const Integer max_itr,
    const std::map<Species, Integer>& max_stoich);

inline std::pair<std::shared_ptr<NetworkModel>, bool> generate_network_from_netfree_model(
    const NetfreeModel& nfm, const std::vector<Species>& seeds, const Integer max_itr)
{
    const std::map<Species, Integer> max_stoich;
    return generate_network_from_netfree_model(
        nfm, seeds, max_itr, max_stoich);
}

} // extras

} // ecell4

#endif /* ECELL4_NETFREE_MODEL_HPP */
