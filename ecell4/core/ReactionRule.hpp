#ifndef ECELL4_REACTION_RULE_HPP
#define ECELL4_REACTION_RULE_HPP

#include "ReactionRuleDescriptor.hpp"

#include "types.hpp"
#include "Quantity.hpp"
#include "Species.hpp"
#include "Attribute.hpp"
#include <stdexcept>
#include <memory>

namespace ecell4
{

class ReactionRule
{
public:

    /**
     * a type of the container of reactants
     * std::multiset allows multiple keys with equal values,
     * but looses the original order at the registration.
     * when changing this type into the ordered one,
     * please modify NetworkModel too.
     */
    typedef std::vector<Species> reactant_container_type;
    typedef std::vector<Species> product_container_type;

public:

    enum policy_type : long
    {
        POLICY_STRICT = 1L << 0,
        POLICY_IMPLICIT = 1L << 1,
        POLICY_DESTROY = 1L << 2
    };

    typedef Attribute::mapped_type attribute_type;

public:

    ReactionRule();
    ReactionRule(const reactant_container_type& reactants,
                 const product_container_type& products);
    ReactionRule(const reactant_container_type& reactants,
                 const product_container_type& products,
                 const Real& k);
    ReactionRule(const reactant_container_type& reactants,
                 const product_container_type& products,
                 const Quantity<Real>& k);
    ReactionRule(const ReactionRule& rr);

    Real k() const;
    void set_k(const Real& k);
    void set_k(const Quantity<Real>& k);
    Quantity<Real> get_k() const;

    const reactant_container_type& reactants() const;
    const product_container_type& products() const;
    void add_reactant(const Species& sp);
    void add_product(const Species& sp);

    const policy_type policy() const;
    void set_policy(const policy_type policy);

    const std::string as_string() const;

    Integer count(const reactant_container_type& reactants) const;
    std::vector<ReactionRule> generate(const reactant_container_type& reactants) const;

    bool has_descriptor() const;
    void set_descriptor(const std::shared_ptr<ReactionRuleDescriptor>& descriptor);
    const std::shared_ptr<ReactionRuleDescriptor>& get_descriptor() const;
    void reset_descriptor();

    /**
     * Attribute
     */
    const Attribute& attributes() const;
    void set_attributes(const Attribute& attr);

    std::vector<std::pair<std::string, attribute_type> > list_attributes() const;
    attribute_type get_attribute(const std::string& key) const;

    template <typename T_>
    T_ get_attribute_as(const std::string& key) const
    {
        return attributes_.get_as<T_>(key);
    }

    template <typename T_>
    void set_attribute(const std::string& key, T_ value)
    {
        attributes_.set<T_>(key, value);
    }

    void remove_attribute(const std::string& key);
    bool has_attribute(const std::string& key) const;

protected:

    Quantity<Real> k_;
    reactant_container_type reactants_;
    product_container_type products_;

    policy_type policy_;
    Attribute attributes_;

    std::shared_ptr<ReactionRuleDescriptor> rr_descriptor_;
    // std::weak_ptr<ReactionRuleDescriptor> rr_descriptor_;
};

inline bool operator<(const ReactionRule& lhs, const ReactionRule& rhs)
{
    if (lhs.reactants() < rhs.reactants())
    {
        return true;
    }
    else if (lhs.reactants() > rhs.reactants())
    {
        return false;
    }
    return (lhs.products() < rhs.products());
}

inline bool operator==(const ReactionRule& lhs, const ReactionRule& rhs)
{
    return ((lhs.reactants() == rhs.reactants())
            && (lhs.products() == rhs.products()));
}

inline bool operator!=(const ReactionRule& lhs, const ReactionRule& rhs)
{
    return !(lhs == rhs);
}

constexpr ReactionRule::policy_type operator|(ReactionRule::policy_type lhs, ReactionRule::policy_type rhs)
{
    return static_cast<ReactionRule::policy_type>(static_cast<long>(lhs) | static_cast<long>(rhs));
}

ReactionRule format_reaction_rule_with_nosort(const ReactionRule& rr);
ReactionRule format_reaction_rule(const ReactionRule& rr);

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

} // ecell4

#endif /* ECELL4_REACTION_RULE_HPP */
