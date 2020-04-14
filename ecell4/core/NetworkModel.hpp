#ifndef ECELL4_NETWORK_MODEL_HPP
#define ECELL4_NETWORK_MODEL_HPP

#include <map>
#include <set>
#include <algorithm>
#include <iterator>
#include <boost/shared_ptr.hpp>

#include "types.hpp"
#include "Species.hpp"
#include "ReactionRule.hpp"
#include "Model.hpp"

#include "Context.hpp"


namespace ecell4
{

class NetworkModel
    : public Model
{
public:

    typedef Model base_type;
    typedef base_type::species_container_type species_container_type;
    typedef base_type::reaction_rule_container_type reaction_rule_container_type;

protected:

    typedef std::map<Species::serial_type,
                     std::vector<reaction_rule_container_type::size_type> >
        first_order_reaction_rules_map_type;
    typedef std::map<std::pair<Species::serial_type, Species::serial_type>,
                     std::vector<reaction_rule_container_type::size_type> >
        second_order_reaction_rules_map_type;

public:

    NetworkModel()
        : base_type(), species_attributes_(), species_attributes_proceed_(), reaction_rules_(),
        first_order_reaction_rules_map_(), second_order_reaction_rules_map_()
    {
        ;
    }

    virtual ~NetworkModel()
    {
        ;
    }

    // ModelTraits

    std::vector<ReactionRule> query_reaction_rules(const Species& sp) const;
    std::vector<ReactionRule> query_reaction_rules(
        const Species& sp1, const Species& sp2) const;

    Integer apply(const Species& pttrn, const Species& sp) const;
    std::vector<ReactionRule> apply(
        const ReactionRule& rr,
        const ReactionRule::reactant_container_type& reactants) const;

    Species apply_species_attributes(const Species& sp) const;

    std::shared_ptr<Model> expand(
        const std::vector<Species>& sp, const Integer max_itr,
        const std::map<Species, Integer>& max_stoich) const
    {
        return std::shared_ptr<Model>(new NetworkModel(*this));
    }

    std::shared_ptr<Model> expand(
        const std::vector<Species>& sp, const Integer max_itr) const
    {
        return std::shared_ptr<Model>(new NetworkModel(*this));
    }

    std::shared_ptr<Model> expand(const std::vector<Species>& sp) const
    {
        return std::shared_ptr<Model>(new NetworkModel(*this));
    }

    // NetworkModelTraits

    bool is_static() const
    {
        return true;
    }

    bool update_species_attribute(const Species& sp);
    void add_species_attribute(const Species& sp, const bool proceed = false);
    bool has_species_attribute(const Species& sp) const;
    void remove_species_attribute(const Species& sp);

    void add_reaction_rule(const ReactionRule& rr);
    void remove_reaction_rule(const ReactionRule& rr);
    bool has_reaction_rule(const ReactionRule& rr) const;

    const reaction_rule_container_type& reaction_rules() const
    {
        return reaction_rules_;
    }

    const species_container_type& species_attributes() const
    {
        return species_attributes_;
    }

    const std::vector<bool>& species_attributes_proceed() const
    {
        return species_attributes_proceed_;
    }

protected:

    void remove_reaction_rule(const reaction_rule_container_type::iterator i);

protected:

    species_container_type species_attributes_;
    std::vector<bool> species_attributes_proceed_;  //XXX:
    reaction_rule_container_type reaction_rules_;

    first_order_reaction_rules_map_type first_order_reaction_rules_map_;
    second_order_reaction_rules_map_type second_order_reaction_rules_map_;
};

} // ecell4

#endif /* ECELL4_NETWORK_MODEL_HPP */
