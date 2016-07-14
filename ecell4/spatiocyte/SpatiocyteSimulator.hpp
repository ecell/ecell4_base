#ifndef __ECELL4_LATTICE_LATTICE_SIMULATOR_HPP
#define __ECELL4_LATTICE_LATTICE_SIMULATOR_HPP

#include <numeric>
#include <boost/shared_ptr.hpp>
#include <boost/lexical_cast.hpp>

#include <ecell4/core/Model.hpp>
#include <ecell4/core/ReactionRule.hpp>
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

    virtual bool check_reaction() const
    {
        return last_reactions().size() > 0;
    }

    const std::vector<SpatiocyteEvent::reaction_type>& last_reactions() const
    {
        return last_event_->reactions();
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

protected:

    boost::shared_ptr<SpatiocyteEvent> create_step_event(
        const Species& species, const Real& t);
    boost::shared_ptr<SpatiocyteEvent> create_zeroth_order_reaction_event(
        const ReactionRule& reaction_rule, const Real& t);
    boost::shared_ptr<SpatiocyteEvent> create_first_order_reaction_event(
        const ReactionRule& reaction_rule, const Real& t);

    void step_();
    void register_events(const Species& species);
    void update_alpha_map();

    void set_last_event_(boost::shared_ptr<const SpatiocyteEvent> event)
    {
        last_event_ = event;
    }

protected:

    scheduler_type scheduler_; boost::shared_ptr<const SpatiocyteEvent> last_event_;
    alpha_map_type alpha_map_;

    Real dt_;
    Real alpha_;
};

} // spatiocyte

} // ecell4

#endif /* __ECELL4_LATTICE_LATTICE_SIMULATOR_HPP */
