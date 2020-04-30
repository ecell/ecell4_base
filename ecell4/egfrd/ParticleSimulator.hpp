#ifndef PARTICLE_SIMULATOR_HPP
#define PARTICLE_SIMULATOR_HPP

#include <boost/bind.hpp>
// #include "Sphere.hpp"
// #include "Cylinder.hpp"
// #include "Box.hpp"
#include "utils/range_support.hpp"
//#include "NetworkRulesWrapper.hpp"
#include "ReactionRuleInfo.hpp"
#include "ReactionRecorder.hpp"
#include "ReactionRecord.hpp"
#include "VolumeClearer.hpp"

#include "NetworkRulesAdapter.hpp"
#include "ReactionRecorderWrapper.hpp"
#include <ecell4/core/SimulatorBase.hpp>
#include <memory>

namespace ecell4
{
namespace egfrd
{

template<typename Tworld_>
struct ParticleSimulatorTraitsBase
{
    typedef Tworld_ world_type;
    typedef Real rate_type;
    typedef Real time_type;
    // typedef int reaction_rule_id_type;
    typedef ecell4::ReactionRule reaction_rule_id_type;
    typedef ReactionRuleInfo<
            reaction_rule_id_type,
            typename world_type::traits_type::species_id_type,
            rate_type> reaction_rule_type;
    //typedef NetworkRulesWrapper<NetworkRules,
    //                           reaction_rule_type> network_rules_type;
    typedef NetworkRulesAdapter<reaction_rule_type> network_rules_type;
    // typedef ReactionRecord<typename world_type::particle_id_type,
    //                        reaction_rule_id_type> reaction_record_type;
    typedef ReactionRecord<typename world_type::particle_id_pair,
                           reaction_rule_id_type> reaction_record_type;
    typedef ReactionRecorder<reaction_record_type> reaction_recorder_type;
    typedef VolumeClearer<typename world_type::particle_shape_type, typename world_type::particle_id_type> volume_clearer_type;

    static const Real minimal_separation_factor();
    static const Real MINIMAL_SEPARATION_FACTOR;
};

template<typename Tworld_>
const Real ParticleSimulatorTraitsBase<Tworld_>::minimal_separation_factor()
{
    return 1 + 1e-7;
}

template<typename Tworld_>
const Real ParticleSimulatorTraitsBase<Tworld_>::MINIMAL_SEPARATION_FACTOR = ParticleSimulatorTraitsBase<Tworld_>::minimal_separation_factor();

template<typename Ttraits_>
class ParticleSimulator;

template<typename Ttraits_>
struct ImmutativeStructureVisitor
{
    typedef Ttraits_ traits_type;
    // typedef typename traits_type::spherical_surface_type spherical_surface_type;
    // typedef typename traits_type::cylindrical_surface_type cylindrical_surface_type;
    // typedef typename traits_type::planar_surface_type planar_surface_type;
    typedef typename traits_type::cuboidal_region_type cuboidal_region_type;

    virtual ~ImmutativeStructureVisitor() {}

    // virtual void operator()(spherical_surface_type const&) const = 0;
    // virtual void operator()(cylindrical_surface_type const&) const = 0;
    // virtual void operator()(planar_surface_type const&) const = 0;
    virtual void operator()(cuboidal_region_type const&) const = 0;
};

template<typename Ttraits_>
struct MutativeStructureVisitor
{
    typedef Ttraits_ traits_type;
    // typedef typename traits_type::spherical_surface_type spherical_surface_type;
    // typedef typename traits_type::cylindrical_surface_type cylindrical_surface_type;
    // typedef typename traits_type::planar_surface_type planar_surface_type;
    typedef typename traits_type::cuboidal_region_type cuboidal_region_type;

    virtual ~MutativeStructureVisitor() {}

    // virtual void operator()(spherical_surface_type&) const = 0;
    // virtual void operator()(cylindrical_surface_type&) const = 0;
    // virtual void operator()(planar_surface_type&) const = 0;
    virtual void operator()(cuboidal_region_type&) const = 0;
};

template<typename Ttraits_>
class ParticleSimulator
    : public ecell4::SimulatorBase<
        typename Ttraits_::world_type,
        typename Ttraits_::world_type::traits_type::model_type>
{
public:

    typedef Ttraits_ traits_type;
    typedef typename traits_type::world_type world_type;
    // typedef Sphere sphere_type;
    // typedef Cylinder cylinder_type;
    // typedef Box box_type;
    // typedef Plane plane_type;

    typedef typename traits_type::network_rules_type network_rules_type;
    typedef typename traits_type::time_type time_type;
    typedef typename traits_type::reaction_record_type reaction_record_type;
    typedef typename traits_type::reaction_recorder_type reaction_recorder_type;
    typedef typename traits_type::volume_clearer_type volume_clearer_type;

    typedef typename world_type::traits_type::rng_type rng_type;
    typedef typename world_type::traits_type::model_type model_type;

    typedef ecell4::SimulatorBase<world_type, model_type> base_type;

public:

    virtual ~ParticleSimulator() {}

    ParticleSimulator(
        const std::shared_ptr<world_type>& world,
        const std::shared_ptr<model_type>& model)
        : base_type(world, model),
        network_rules_(new network_rules_type(model)),
        rrec_(new ReactionRecorderWrapper<reaction_record_type>()),
        dt_(0.), paranoiac_(false)
    {
        ;
    }

    ParticleSimulator(
        const std::shared_ptr<world_type>& world)
        : base_type(world),
        network_rules_(new network_rules_type(this->model())),
        rrec_(new ReactionRecorderWrapper<reaction_record_type>()),
        dt_(0.), paranoiac_(false)
    {
        ;
    }

    std::shared_ptr<network_rules_type const> const& network_rules() const
    {
        return network_rules_;
    }

    // std::shared_ptr<reaction_recorder_type> const& reaction_recorder() const
    // {
    //     return rrec_;
    // }

    // std::shared_ptr<reaction_recorder_type>& reaction_recorder()
    // {
    //     return rrec_;
    // }

    inline rng_type& rng() const
    {
        return (*(*base_type::world_).rng().get());
    }

    virtual void set_dt(const Real& dt)
    {
        std::cerr << "WARN: set_dt(const Real&) was just ignored." << std::endl;
        dt_ = dt;
    }

    virtual time_type dt() const
    {
        return dt_; //XXX: dt has no mean for egfrd.
    }

    bool const& paranoiac() const
    {
        return paranoiac_;
    }

    bool& paranoiac()
    {
        return paranoiac_;
    }

    void set_paranoiac(const bool val)
    {
        paranoiac_ = val;
    }

    virtual void step() = 0;
    virtual bool step(const time_type& upto) = 0;

protected:
    std::shared_ptr<network_rules_type const> network_rules_;
    std::shared_ptr<reaction_recorder_type> rrec_;
    time_type dt_;
    bool paranoiac_;

};

} // egfrd
} // ecell4
#endif /* PARTICLE_SIMULATOR_HPP */
