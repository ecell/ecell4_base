#ifndef __ECELL4_SPATIOCYTE_EVENT_HPP
#define __ECELL4_SPATIOCYTE_EVENT_HPP

#include <ecell4/core/ReactionRule.hpp>
#include <ecell4/core/EventScheduler.hpp>
#include <ecell4/core/VoxelPool.hpp>

namespace ecell4
{

namespace spatiocyte {

class SpatiocyteSimulator;

class ReactionInfo
{
public:

    typedef std::pair<ParticleID, Voxel> particle_id_pair_type;
    typedef std::vector<particle_id_pair_type> container_type;

public:

    ReactionInfo() : t_(0), reactants_(), products_() {}

    ReactionInfo(const Real t) : t_(t), reactants_(), products_() {}

    ReactionInfo(
        const Real t,
        const container_type& reactants,
        const container_type& products)
        : t_(t), reactants_(reactants), products_(products) {}

    ReactionInfo(const ReactionInfo& another)
        : t_(another.t()), reactants_(another.reactants()), products_(another.products()) {}

    Real t() const
    {
        return t_;
    }

    const container_type& reactants() const
    {
        return reactants_;
    }

    void add_reactant(const particle_id_pair_type& pid_pair)
    {
        reactants_.push_back(pid_pair);
    }

    const container_type& products() const
    {
        return products_;
    }

    void add_product(const particle_id_pair_type& pid_pair)
    {
        products_.push_back(pid_pair);
    }

protected:

    Real t_;
    container_type reactants_, products_;
};

struct SpatiocyteEvent : public Event
{
public:
    typedef ReactionInfo reaction_info_type;
    typedef std::pair<ReactionRule, reaction_info_type> reaction_type;

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
    StepEvent(SpatiocyteSimulator* sim, const Species& species, const Real& t,
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

    void walk(const Real& alpha) const;

protected:

    SpatiocyteSimulator* sim_;
    Species species_;
    VoxelPool* mt_;
    const Real alpha_;
};

struct ZerothOrderReactionEvent : SpatiocyteEvent
{
    ZerothOrderReactionEvent(
        SpatiocyteSimulator* sim, const ReactionRule& rule, const Real& t);

    virtual ~ZerothOrderReactionEvent() {}
    virtual void fire();

    Real draw_dt();
    virtual void interrupt(Real const& t)
    {
        time_ = t + draw_dt();
    }

protected:

    SpatiocyteSimulator* sim_;
    ReactionRule rule_;
};

struct FirstOrderReactionEvent : SpatiocyteEvent
{
    FirstOrderReactionEvent(
        SpatiocyteSimulator* sim, const ReactionRule& rule, const Real& t);

    virtual ~FirstOrderReactionEvent() {}
    virtual void fire();

    Real draw_dt();
    virtual void interrupt(Real const& t)
    {
        time_ = t + draw_dt();
    }

protected:

    SpatiocyteSimulator* sim_;
    ReactionRule rule_;
};

} // spatiocyte

} // ecell4

#endif /* __ECELL4_SPATIOCYTE_EVENT_HPP */
