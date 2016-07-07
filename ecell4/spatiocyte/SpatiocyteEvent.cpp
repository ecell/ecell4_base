#include "SpatiocyteEvent.hpp"
#include "SpatiocyteSimulator.hpp" // TODO: remove this include line

namespace ecell4 {

namespace spatiocyte {

/// StepEvent

StepEvent::StepEvent(SpatiocyteSimulator* sim, const Species& species, const Real& t,
    const Real alpha)
    : SpatiocyteEvent(t), sim_(sim), world_(sim->world()), species_(species), alpha_(alpha)
{
    const SpatiocyteWorld::molecule_info_type
        minfo(world_->get_molecule_info(species));
    const Real R(minfo.radius);
    const Real D(minfo.D);
    const VoxelPool* mtype(world_->find_voxel_pool(species));
    // const Real R(world_->voxel_radius());
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

    nids_.clear();
    for (unsigned int i(0); i < 12; ++i)
        nids_.push_back(i);
}

void StepEvent::fire()
{
    walk(alpha_);
    time_ += dt_;
}

void StepEvent::walk(const Real& alpha)
{
    if (alpha < 0 || alpha > 1)
    {
        return; // INVALID ALPHA VALUE
    }

    const boost::shared_ptr<RandomNumberGenerator>& rng(world_->rng());
    const MoleculePool* mtype(world_->find_molecule_pool(species_));

    if (mtype->get_dimension() == Shape::THREE)
        walk_in_space_(mtype, alpha);
    else // dimension == TWO, etc.
        walk_on_surface_(mtype, alpha);
}

void StepEvent::walk_in_space_(const MoleculePool* mtype, const Real& alpha)
{
    const boost::shared_ptr<RandomNumberGenerator>& rng(world_->rng());
    MoleculePool::container_type voxels;
    copy(mtype->begin(), mtype->end(), back_inserter(voxels));

    std::size_t idx(0);
    for (MoleculePool::container_type::iterator itr(voxels.begin());
         itr != voxels.end(); ++itr)
    {
        const Integer rnd(rng->uniform_int(0, 11));
        const SpatiocyteWorld::coordinate_id_pair_type& info(*itr);
        if (world_->find_voxel_pool(info.coordinate) != mtype)
        {
            // should skip if a voxel is not the target species.
            // when reaction has occured before, a voxel can be changed.
            continue;
        }
        const SpatiocyteWorld::coordinate_type neighbor(
                world_->get_neighbor_boundary(info.coordinate, rnd));
        if (world_->can_move(info.coordinate, neighbor))
        {
            if (rng->uniform(0,1) <= alpha)
                world_->move(info.coordinate, neighbor, /*candidate=*/idx);
        }
        else
        {
            sim_->attempt_reaction_(info, neighbor, alpha);
        }
        ++idx;
    }
}

void StepEvent::walk_on_surface_(const MoleculePool* mtype, const Real& alpha)
{
    const boost::shared_ptr<RandomNumberGenerator>& rng(world_->rng());
    MoleculePool::container_type voxels;
    copy(mtype->begin(), mtype->end(), back_inserter(voxels));

    const VoxelPool* location(mtype->location());
    std::size_t idx(0);
    for (MoleculePool::container_type::iterator itr(voxels.begin());
         itr != voxels.end(); ++itr)
    {
        const SpatiocyteWorld::coordinate_id_pair_type& info(*itr);
        if (world_->find_voxel_pool(info.coordinate) != mtype)
        {
            // should skip if a voxel is not the target species.
            // when reaction has occured before, a voxel can be changed.
            continue;
        }

        ecell4::shuffle(*(rng.get()), nids_);
        for (std::vector<unsigned int>::const_iterator itr(nids_.begin());
             itr != nids_.end(); ++itr)
        {
            const SpatiocyteWorld::coordinate_type neighbor(
                    world_->get_neighbor_boundary(info.coordinate, *itr));
            const VoxelPool* target(world_->find_voxel_pool(neighbor));

            if (target->get_dimension() > mtype->get_dimension())
                continue;

            if (world_->can_move(info.coordinate, neighbor))
            {
                if (rng->uniform(0,1) <= alpha)
                    world_->move(info.coordinate, neighbor, /*candidate=*/idx);
                break;
            }
            else if (sim_->attempt_reaction_(info, neighbor, alpha).first != SpatiocyteSimulator::NO_REACTION)
                break;
        }
        ++idx;
    }
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
