#include "SpatiocyteEvent.hpp"
#include "SpatiocyteSimulator.hpp" // TODO: remove this include line

namespace ecell4 {

namespace spatiocyte {

/// ZerothOrderReactionEvent

ZerothOrderReactionEvent::ZerothOrderReactionEvent(
    SpatiocyteSimulator* sim, const ReactionRule& rule, const Real& t)
    : SpatiocyteEvent(t), sim_(sim), rule_(rule)
{
    time_ = t + draw_dt();
}

void ZerothOrderReactionEvent::fire()
{
    const std::pair<bool, reaction_type> reaction(
            sim_->apply_zeroth_order_reaction_(rule_));
    if (reaction.first)
        push_reaction(reaction.second);
    time_ += draw_dt();
}

Real ZerothOrderReactionEvent::draw_dt()
{
    const Real k(rule_.k());
    const Real p = k * sim_->world()->volume();
    Real dt(inf);
    if (p != 0.)
    {
        const Real rnd(sim_->world()->rng()->uniform(0.,1.));
        dt = - log(1 - rnd) / p;
    }
    return dt;
}


/// FirstOrderReactionEvent

FirstOrderReactionEvent::FirstOrderReactionEvent(
    SpatiocyteSimulator* sim, const ReactionRule& rule, const Real& t)
    : SpatiocyteEvent(t), sim_(sim), rule_(rule)
{
    //assert(rule_.reactants().size() == 1);
    time_ = t + draw_dt();
}

void FirstOrderReactionEvent::fire()
{
    const Species& reactant(*(rule_.reactants().begin()));
    const std::pair<bool, reaction_type> reaction(
            sim_->apply_first_order_reaction_(rule_, sim_->world()->choice(reactant)));
    if (reaction.first)
        push_reaction(reaction.second);
    time_ += draw_dt();
}

Real FirstOrderReactionEvent::draw_dt()
{
    const Species& reactant(*(rule_.reactants().begin()));
    const Integer num_r(sim_->world()->num_voxels_exact(reactant));
    const Real k(rule_.k());
    const Real p = k * num_r;
    Real dt(inf);
    if (p != 0.)
    {
        const Real rnd(sim_->world()->rng()->uniform(0.,1.));
        dt = - log(1 - rnd) / p;
    }
    return dt;
}

} // spatiocyte

} // ecell4
