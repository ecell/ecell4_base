#ifndef __ECELL4_LATTICE_LATTICE_SIMULATOR_HPP
#define __ECELL4_LATTICE_LATTICE_SIMULATOR_HPP

#include <numeric>
#include <boost/shared_ptr.hpp>
#include <boost/lexical_cast.hpp>

#include <ecell4/core/Model.hpp>
#include <ecell4/core/ReactionRule.hpp>
// #include <ecell4/core/Reaction.hpp>
#include <ecell4/core/VoxelPool.hpp>
#include <ecell4/core/SimulatorBase.hpp>
#include <ecell4/core/RandomNumberGenerator.hpp>
#include <ecell4/core/EventScheduler.hpp>
#include <ecell4/core/get_mapper_mf.hpp>

#include "SpatiocyteWorld.hpp"
#include "SpatiocyteEvent.hpp"

namespace ecell4
{

namespace spatiocyte
{

class SpatiocyteSimulator
    : public SimulatorBase<Model, SpatiocyteWorld>
{
public:

    typedef SimulatorBase<Model, SpatiocyteWorld> base_type;

    typedef SpatiocyteEvent::reaction_type reaction_type;

    typedef EventSchedulerBase<SpatiocyteEvent> scheduler_type;

    typedef utils::get_mapper_mf<Species, Real>::type alpha_map_type;

public:

    SpatiocyteSimulator(
            boost::shared_ptr<Model> model,
            boost::shared_ptr<SpatiocyteWorld> world,
            const Real alpha = 1.0)
        : base_type(model, world), alpha_(alpha)
    {
        initialize();
    }

    SpatiocyteSimulator(
            boost::shared_ptr<SpatiocyteWorld> world,
            const Real alpha = 1.0)
        : base_type(world), alpha_(alpha)
    {
        initialize();
    }

    virtual Real dt() const
    {
        return dt_;
    }

    void initialize();
    void finalize();
    void step();
    bool step(const Real& upto);
    Real calculate_alpha(const ReactionRule& rule) const;

    virtual bool check_reaction() const
    {
        return last_reactions_.size() > 0;
    }

    std::vector<std::pair<ReactionRule, ReactionInfo> > last_reactions() const
    {
        return last_reactions_;
    }

    void set_alpha(const Real alpha)
    {
        if (alpha < 0 || alpha > 1)
        {
            return;  // XXX: ValueError
        }

        alpha_ = alpha;
        initialize();
    }

    Real get_alpha() const
    {
        return alpha_;
    }

    // TODO: remove the below public functions
    static Real calculate_dimensional_factor(const VoxelPool* mt0, const VoxelPool* mt1,
            boost::shared_ptr<SpatiocyteWorld> world);

protected:

    boost::shared_ptr<SpatiocyteEvent> create_step_event(
        const Species& species, const Real& t);
    boost::shared_ptr<SpatiocyteEvent> create_zeroth_order_reaction_event(
        const ReactionRule& reaction_rule, const Real& t);
    boost::shared_ptr<SpatiocyteEvent> create_first_order_reaction_event(
        const ReactionRule& reaction_rule, const Real& t);

    // void register_reactant_species(
    //     const SpatiocyteWorld::coordinate_id_pair_type pinfo, reaction_type& reaction) const;

    void step_();
    void register_events(const Species& species);
    // void register_step_event(const Species& species);
    void update_alpha_map();

    const std::string get_serial(
        const SpatiocyteWorld::coordinate_type coord) const
    {
        const VoxelPool* mtype(world_->find_voxel_pool(coord));
        return mtype->is_vacant() ? "" : mtype->species().serial();
    }

    const std::string get_location(
        const SpatiocyteWorld::coordinate_type coord) const
    {
        const VoxelPool* mtype(world_->find_voxel_pool(coord));
        if (mtype->is_vacant())
        {
            return "";
        }
        const VoxelPool* ltype(mtype->location());
        return ltype->is_vacant() ? "" : ltype->species().serial();
    }

protected:

    scheduler_type scheduler_;
    std::vector<reaction_type> last_reactions_;
    alpha_map_type alpha_map_;
    //map<Species> alpha_map_;

    Real dt_;
    Real alpha_;
};

} // spatiocyte

} // ecell4

#endif /* __ECELL4_LATTICE_LATTICE_SIMULATOR_HPP */
