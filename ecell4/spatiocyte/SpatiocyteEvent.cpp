#include "SpatiocyteEvent.hpp"
#include "SpatiocyteSimulator.hpp" // TODO: remove this include line

namespace ecell4 {

namespace spatiocyte {

/// StepEvent

StepEvent::StepEvent(SpatiocyteSimulator* sim, const Species& species, const Real& t,
    const Real alpha)
    : SpatiocyteEvent(t), sim_(sim), species_(species), alpha_(alpha)
{
    const SpatiocyteWorld::molecule_info_type
        minfo(sim_->world()->get_molecule_info(species));
    const Real R(minfo.radius);
    const Real D(minfo.D);
    const VoxelPool* mtype(sim_->world()->find_voxel_pool(species));
    // const Real R(sim_->world()->voxel_radius());
    // Real D = boost::lexical_cast<Real>(species.get_attribute("D"));
    if (D <= 0)
    {
        dt_ = inf;
    } else if(mtype->get_dimension() == Shape::THREE) {
        dt_ = 2 * R * R / 3 / D * alpha_;
    } else if(mtype->get_dimension() == Shape::TWO) {
        // TODO: Regular Lattice
        // dt_  = pow((2*sqrt(2.0)+4*sqrt(3.0)+3*sqrt(6.0)+sqrt(22.0))/
        //           (6*sqrt(2.0)+4*sqrt(3.0)+3*sqrt(6.0)), 2) * R * R / D * alpha_;
        dt_ = R * R / D * alpha_;
    } else if(mtype->get_dimension() == Shape::ONE) {
        dt_ = 2 * R * R / D * alpha_;
    }
    else
    {
        throw NotSupported(
            "The dimension of a structure must be two or three.");
    }

    time_ = t + dt_;
    // time_ = t;
}

void StepEvent::fire()
{
    sim_->walk(species_, alpha_);
    time_ += dt_;
}


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
