#ifndef PARTICLE_SIMULATOR_HPP
#define PARTICLE_SIMULATOR_HPP

#include <boost/shared_ptr.hpp>
#include "Sphere.hpp"
#include "Cylinder.hpp"
#include "Box.hpp"
#include "CuboidalRegion.hpp"
#include "PlanarSurface.hpp"
#include "CylindricalSurface.hpp"
#include "SphericalSurface.hpp"
#include "NetworkRules.hpp"
#include "NetworkRulesWrapper.hpp"
#include "ReactionRuleInfo.hpp"
#include "ReactionRecorder.hpp"
#include "ReactionRecord.hpp"

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
    typedef Plane<typename world_type::length_type> plane_type;
    typedef SphericalSurface<typename world_type::traits_type> spherical_surface_type;
    typedef CylindricalSurface<typename world_type::traits_type> cylindrical_surface_type;
    typedef PlanarSurface<typename world_type::traits_type> planar_surface_type;
    typedef CuboidalRegion<typename world_type::traits_type> cuboidal_region_type;
    typedef ReactionRuleInfo<
            reaction_rule_id_type,
            typename world_type::traits_type::species_id_type,
            rate_type> reaction_rule_type;
    typedef NetworkRulesWrapper<NetworkRules,
                                reaction_rule_type> network_rules_type;
    typedef ReactionRecord<typename world_type::particle_id_type,
                           reaction_rule_id_type> reaction_record_type;
    typedef ReactionRecorder<reaction_record_type> reaction_recorder_type;
};

template<typename Ttraits_>
class ParticleSimulator
{
public:
    typedef Ttraits_ traits_type;
    typedef typename traits_type::world_type world_type;
    typedef typename traits_type::network_rules_type network_rules_type;
    typedef typename world_type::traits_type::rng_type rng_type;
    typedef typename traits_type::time_type time_type;
    typedef typename traits_type::reaction_record_type reaction_record_type;
    typedef typename traits_type::reaction_recorder_type reaction_recorder_type;

public:
    virtual ~ParticleSimulator() {}

    ParticleSimulator(boost::shared_ptr<world_type> world,
                      boost::shared_ptr<network_rules_type const> network_rules,
                      rng_type& rng, reaction_recorder_type& rrec)
        : world_(world), network_rules_(network_rules), rng_(rng), rrec_(rrec),
          t_(0.), dt_(0.), num_steps_(0) {}

    boost::shared_ptr<world_type> world() const
    {
        return world_;
    }

    boost::shared_ptr<network_rules_type const> network_rules() const
    {
        return network_rules_;
    }

    rng_type& rng() const
    {
        return rng_;
    }

    time_type t() const
    {
        return t_;
    }

    time_type dt() const
    {
        return dt_;
    }

    int num_steps() const
    {
        return num_steps_;
    }

    virtual void step() = 0;

    virtual bool step(time_type const& upto) = 0;

protected:
    boost::shared_ptr<world_type> world_;
    boost::shared_ptr<network_rules_type const> network_rules_;
    rng_type& rng_;
    reaction_recorder_type& rrec_;
    time_type t_;
    time_type dt_;
    int num_steps_;
};

#endif /* PARTICLE_SIMULATOR_HPP */
