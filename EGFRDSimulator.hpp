#ifndef EGFRDSIMULATOR_HPP
#define EGFRDSIMULATOR_HPP

#include "ShellID.hpp"
#include "DomainID.hpp"
#include "Shell.hpp"
#include "Cylinder.hpp"
#include "Sphere.hpp"
#include "NetworkRulesWrapper.hpp"
#include "ReactionRuleInfo.hpp"
#include "GSLRandomNumberGenerator.hpp"

template<typename Tworld_, typename Tmodel_>
struct EGFRDSimulatorTraitsBase
{
    typedef Tworld_ world_type;
    typedef Tmodel_ model_type;
    typedef ShellID shell_id_type;
    typedef DomainID domain_id_type;
    typedef Real rate_type;
    typedef Real time_type;
    typedef int reaction_rule_id_type;
    typedef Sphere<typename world_type::length_type> sphere_type;
    typedef Cylinder<typename world_type::length_type> cylinder_type;
    typedef Box<typename world_type::length_type> box_type;
    typedef Shell<sphere_type, domain_id_type> spherical_shell_type;
    typedef Shell<cylinder_type, domain_id_type> cylindrical_shell_type;
    typedef ReactionRuleInfo<
            reaction_rule_id_type,
            typename world_type::traits_type::species_id_type,
            rate_type> reaction_rule_type;
    typedef NetworkRulesWrapper<typename model_type::network_rules_type,
                                reaction_rule_type> network_rules_type;
    typedef GSLRandomNumberGenerator rng_type;
};

#endif /* EGFRDSIMULATOR_HPP */
