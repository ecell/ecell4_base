#ifndef ECELL4_LATTICE_LATTICE_SIMULATOR_HPP
#define ECELL4_LATTICE_LATTICE_SIMULATOR_HPP

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
            boost::shared_ptr<SpatiocyteWorld> world)
        : base_type(model, world)
    {
        initialize();
    }

    SpatiocyteSimulator(
            boost::shared_ptr<SpatiocyteWorld> world)
        : base_type(world)
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
        // return last_event_->reactions();
        return last_reactions_;
    }

protected:

    boost::shared_ptr<SpatiocyteEvent> create_step_event(
        const Species& species, const Real& t, const Real& alpha);
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

    std::vector<reaction_type> last_reactions_;

    Real dt_;
};

} // spatiocyte

} // ecell4

#endif /* ECELL4_LATTICE_LATTICE_SIMULATOR_HPP */
