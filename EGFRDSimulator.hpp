#ifndef EGFRDSIMULATOR_HPP
#define EGFRDSIMULATOR_HPP

#include "ShellID.hpp"
#include "DomainID.hpp"
#include "SphericalShell.hpp"
#include "CylindricalShell.hpp"
#include "NetworkRulesWrapper.hpp"
#include "ReactionRuleInfo.hpp"

template<typename Tworld_, typename Tmodel_>
struct EGFRDSimulatorTraitsBase
{
    typedef Tworld_ world_type;
    typedef Tmodel_ model_type;
    typedef ShellID shell_id_type;
    typedef DomainID domain_id_type;
    typedef Real rate_type;
    typedef int reaction_rule_id_type;
    typedef SphericalShell<typename world_type::length_type,
                           domain_id_type> spherical_shell_type;
    typedef CylindricalShell<typename world_type::length_type,
                             domain_id_type> cylindrical_shell_type;
    typedef ReactionRuleInfo<
            reaction_rule_id_type,
            typename world_type::traits_type::species_id_type,
            rate_type> reaction_rule_type;
    typedef NetworkRulesWrapper<typename model_type::network_rules_type,
                                reaction_rule_type> network_rules_type;
};

#endif /* EGFRDSIMULATOR_HPP */
