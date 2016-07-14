#ifndef __ECELL4_SPATIOCYTE_EVENT_HPP
#define __ECELL4_SPATIOCYTE_EVENT_HPP

#include <ecell4/core/ReactionRule.hpp>
#include <ecell4/core/EventScheduler.hpp>
#include <ecell4/core/Model.hpp>
#include "SpatiocyteReactions.hpp"
#include "SpatiocyteWorld.hpp"

namespace ecell4
{

namespace spatiocyte
{

class SpatiocyteSimulator;

struct SpatiocyteEvent : public Event
{
public:
    typedef std::pair<ReactionRule, ReactionInfo> reaction_type;

    SpatiocyteEvent(Real const& time) : Event(time) {}
    virtual ~SpatiocyteEvent() {}

    std::vector<reaction_type> last_reactions() const
    {
        return last_reactions_;
    }

    void push_reaction(const reaction_type& reaction)
    {
        last_reactions_.push_back(reaction);
    }

protected:
    std::vector<reaction_type> last_reactions_;

};

struct StepEvent : SpatiocyteEvent
{
    StepEvent(SpatiocyteSimulator* sim, boost::shared_ptr<Model> model,
            boost::shared_ptr<SpatiocyteWorld> world, const Species& species, const Real& t,
        const Real alpha=1.0);

    virtual ~StepEvent() {}
    virtual void fire();

    Species const& species() const
    {
        return species_;
    }

    Real const& alpha() const
    {
        return alpha_;
    }

    void walk(const Real& alpha);

protected:

    typedef enum
    {
        NO_REACTION = 0,
        REACTION_FAILED = 1,
        REACTION_SUCCEEDED = 2
    } attempt_reaction_result_type;

    void walk_in_space_(const MoleculePool* mtype, const Real& alpha);
    void walk_on_surface_(const MoleculePool* mtype, const Real& alpha);
    std::pair<attempt_reaction_result_type, reaction_type> attempt_reaction_(
        const SpatiocyteWorld::coordinate_id_pair_type& info,
        const SpatiocyteWorld::coordinate_type to_coord, const Real& alpha);

    SpatiocyteSimulator* sim_;
    boost::shared_ptr<Model> model_;
    boost::shared_ptr<SpatiocyteWorld> world_;
    Species species_;
    VoxelPool* mt_;
    const Real alpha_;
    std::vector<unsigned int> nids_; // neighbor indexes
};

struct ZerothOrderReactionEvent : SpatiocyteEvent
{
    ZerothOrderReactionEvent(
        boost::shared_ptr<SpatiocyteWorld> world, const ReactionRule& rule, const Real& t);

    virtual ~ZerothOrderReactionEvent() {}
    virtual void fire();

    Real draw_dt();
    virtual void interrupt(Real const& t)
    {
        time_ = t + draw_dt();
    }

protected:

    boost::shared_ptr<SpatiocyteWorld> world_;
    ReactionRule rule_;
};

struct FirstOrderReactionEvent : SpatiocyteEvent
{
    FirstOrderReactionEvent(
        boost::shared_ptr<SpatiocyteWorld> world, const ReactionRule& rule, const Real& t);

    virtual ~FirstOrderReactionEvent() {}
    virtual void fire();

    Real draw_dt();
    virtual void interrupt(Real const& t)
    {
        time_ = t + draw_dt();
    }

protected:

    boost::shared_ptr<SpatiocyteWorld> world_;
    ReactionRule rule_;
};

} // spatiocyte

} // ecell4

#endif /* __ECELL4_SPATIOCYTE_EVENT_HPP */
