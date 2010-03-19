#ifndef PARTICLE_SIMULATOR_HPP
#define PARTICLE_SIMULATOR_HPP

#include "Sphere.hpp"
#include "Cylinder.hpp"
#include "Box.hpp"
#include "Surface.hpp"
#include "Region.hpp"
#include "NetworkRules.hpp"
#include "NetworkRulesWrapper.hpp"
#include "ReactionRuleInfo.hpp"
#include "GSLRandomNumberGenerator.hpp"

template<typename Tworld_>
struct ParticleSimulatorTraitsBase
{
    typedef Tworld_ world_type;
    typedef Real rate_type;
    typedef Real time_type;
    typedef int reaction_rule_id_type;
    typedef Sphere<typename world_type::length_type> sphere_type;
    typedef Cylinder<typename world_type::length_type> cylinder_type;
    typedef Box<typename world_type::length_type> box_type;
    typedef Surface<typename world_type::surface_id_type,
                    sphere_type> spherical_surface_type;
    typedef Surface<typename world_type::surface_id_type,
                    cylinder_type> cylindrical_surface_type;
    typedef Surface<typename world_type::surface_id_type,
                    box_type> planar_surface_type;
    typedef Region<typename world_type::surface_id_type,
                    box_type> cuboidal_region_type;
    typedef ReactionRuleInfo<
            reaction_rule_id_type,
            typename world_type::traits_type::species_id_type,
            rate_type> reaction_rule_type;
    typedef NetworkRulesWrapper<NetworkRules,
                                reaction_rule_type> network_rules_type;
    typedef GSLRandomNumberGenerator rng_type;
};

#endif /* PARTICLE_SIMULATOR_HPP */
