#ifndef PARTICLE_SIMULATOR_HPP
#define PARTICLE_SIMULATOR_HPP

#include <boost/bind.hpp>
#include <boost/shared_ptr.hpp>
#include "Sphere.hpp"
#include "Cylinder.hpp"
#include "Box.hpp"
#include "ParticleSimulationStructure.hpp"
#include "CuboidalRegion.hpp"
#include "PlanarSurface.hpp"
#include "CylindricalSurface.hpp"
#include "SphericalSurface.hpp"
#include "NetworkRules.hpp"
#include "NetworkRulesWrapper.hpp"
#include "ReactionRuleInfo.hpp"
#include "ReactionRecorder.hpp"
#include "ReactionRecord.hpp"
#include "VolumeClearer.hpp"

template<typename Tworld_>
struct ParticleSimulatorTraitsBase
{
    typedef Tworld_ world_type;
    typedef Real rate_type;
    typedef Real time_type;
    typedef int reaction_rule_id_type;
    typedef ReactionRuleInfo<
            reaction_rule_id_type,
            typename world_type::traits_type::species_id_type,
            rate_type> reaction_rule_type;
    typedef NetworkRulesWrapper<NetworkRules,
                                reaction_rule_type> network_rules_type;
    typedef ReactionRecord<typename world_type::particle_id_type,
                           reaction_rule_id_type> reaction_record_type;
    typedef ReactionRecorder<reaction_record_type> reaction_recorder_type;
    typedef VolumeClearer<typename world_type::particle_shape_type, typename world_type::particle_id_type> volume_clearer_type;

    static const Real MINIMAL_SEPARATION_FACTOR = (1.0 + 1e-7);
};

template<typename Ttraits_>
class ParticleSimulator;

template<typename Ttraits_>
struct ImmutativeStructureVisitor
{
    typedef Ttraits_ traits_type;
    typedef typename ParticleSimulator<traits_type>::spherical_surface_type spherical_surface_type;
    typedef typename ParticleSimulator<traits_type>::cylindrical_surface_type cylindrical_surface_type;
    typedef typename ParticleSimulator<traits_type>::planar_surface_type planar_surface_type;
    typedef typename ParticleSimulator<traits_type>::cuboidal_region_type cuboidal_region_type;

    virtual ~ImmutativeStructureVisitor() {}

    virtual void operator()(spherical_surface_type const&) const = 0;

    virtual void operator()(cylindrical_surface_type const&) const = 0;

    virtual void operator()(planar_surface_type const&) const = 0;

    virtual void operator()(cuboidal_region_type const&) const = 0;
};

template<typename Ttraits_>
struct MutativeStructureVisitor
{
    typedef Ttraits_ traits_type;
    typedef typename ParticleSimulator<traits_type>::spherical_surface_type spherical_surface_type;
    typedef typename ParticleSimulator<traits_type>::cylindrical_surface_type cylindrical_surface_type;
    typedef typename ParticleSimulator<traits_type>::planar_surface_type planar_surface_type;
    typedef typename ParticleSimulator<traits_type>::cuboidal_region_type cuboidal_region_type;

    virtual ~MutativeStructureVisitor() {}

    virtual void operator()(spherical_surface_type&) const = 0;

    virtual void operator()(cylindrical_surface_type&) const = 0;

    virtual void operator()(planar_surface_type&) const = 0;

    virtual void operator()(cuboidal_region_type&) const = 0;
};

template<typename Ttraits_>
class ParticleSimulator
{
public:
    typedef Ttraits_ traits_type;
    typedef typename traits_type::world_type world_type;
    typedef Sphere<typename world_type::length_type> sphere_type;
    typedef Cylinder<typename world_type::length_type> cylinder_type;
    typedef Box<typename world_type::length_type> box_type;
    typedef Plane<typename world_type::length_type> plane_type;
    typedef ParticleSimulationStructure<traits_type> particle_simulation_structure_type;
    typedef Surface<traits_type> surface_type;
    typedef Region<traits_type> region_type;
    typedef SphericalSurface<traits_type> spherical_surface_type;
    typedef CylindricalSurface<traits_type> cylindrical_surface_type;
    typedef PlanarSurface<traits_type> planar_surface_type;
    typedef CuboidalRegion<traits_type> cuboidal_region_type;
    typedef typename traits_type::network_rules_type network_rules_type;
    typedef typename world_type::traits_type::rng_type rng_type;
    typedef typename traits_type::time_type time_type;
    typedef typename traits_type::reaction_record_type reaction_record_type;
    typedef typename traits_type::reaction_recorder_type reaction_recorder_type;
    typedef typename traits_type::volume_clearer_type volume_clearer_type;

public:
    virtual ~ParticleSimulator() {}

    ParticleSimulator(boost::shared_ptr<world_type> world,
                      boost::shared_ptr<network_rules_type const> network_rules,
                      rng_type& rng)
        : world_(world), network_rules_(network_rules), rrec_(), rng_(rng),
          t_(0.), dt_(0.), num_steps_(0), paranoiac_(false) {}

    boost::shared_ptr<world_type> const& world() const
    {
        return world_;
    }

    boost::shared_ptr<network_rules_type const> const& network_rules() const
    {
        return network_rules_;
    }

    boost::shared_ptr<reaction_recorder_type> const& reaction_recorder() const
    {
        return rrec_;
    }

    boost::shared_ptr<reaction_recorder_type>& reaction_recorder()
    {
        return rrec_;
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

    bool const& paranoiac() const
    {
        return paranoiac_;
    }

    bool& paranoiac()
    {
        return paranoiac_;
    }

    int num_steps() const
    {
        return num_steps_;
    }

    virtual void step() = 0;

    virtual bool step(time_type upto) = 0;

protected:
    boost::shared_ptr<world_type> world_;
    boost::shared_ptr<network_rules_type const> network_rules_;
    boost::shared_ptr<reaction_recorder_type> rrec_;
    rng_type& rng_;
    time_type t_;
    time_type dt_;
    int num_steps_;
    bool paranoiac_;
};

#endif /* PARTICLE_SIMULATOR_HPP */
