#include "SpatiocyteEvent.hpp"

namespace ecell4 {

namespace spatiocyte {

/// ZerothOrderReactionEvent

ZerothOrderReactionEvent::ZerothOrderReactionEvent(
    boost::shared_ptr<SpatiocyteWorld> world, const ReactionRule& rule, const Real& t)
    : SpatiocyteEvent(t), world_(world), rule_(rule)
{
    time_ = t + draw_dt();
}

void ZerothOrderReactionEvent::fire_()
{
    ReactionInfo rinfo(world_->t());

    for (ReactionRule::product_container_type::const_iterator
        i(rule_.products().begin());
        i != rule_.products().end(); ++i)
    {
        const Species& sp(*i);
        const SpatiocyteWorld::molecule_info_type info(world_->get_molecule_info(sp));

        if (boost::shared_ptr<VoxelPool> location = world_->find_voxel_pool(Species(info.loc)))
        {
            if (location->size() == 0)
            {
                time_ =+ draw_dt();
                return;
            }

            while (true)
            {
                const SpatiocyteWorld::coordinate_type
                    coord(world_->rng()->uniform_int(0, world_->size() - 1));

                if (world_->get_voxel_pool_at(coord) != location)
                {
                    continue;
                }

                if (boost::optional<ParticleID> pid = world_->new_voxel(sp, coord))
                {
                    rinfo.add_product(ReactionInfo::Item(*pid, sp, coord));
                    break;
                }
            }
        }
    }
    push_reaction(std::make_pair(rule_, rinfo));
    time_ += draw_dt();
}

Real ZerothOrderReactionEvent::draw_dt()
{
    const Real k(rule_.k());
    //XXX: const Real p = k * world_->volume();
    const Real p = k * world_->actual_volume();
    Real dt(inf);
    if (p != 0.)
    {
        const Real rnd(world_->rng()->uniform(0.,1.));
        dt = - log(1 - rnd) / p;
    }
    return dt;
}


/// FirstOrderReactionEvent

FirstOrderReactionEvent::FirstOrderReactionEvent(
    boost::shared_ptr<SpatiocyteWorld> world, const ReactionRule& rule, const Real& t)
    : SpatiocyteEvent(t), world_(world), rule_(rule)
{
    //assert(rule_.reactants().size() == 1);
    time_ = t + draw_dt();
}

void FirstOrderReactionEvent::fire_()
{
    const ReactionInfo::Item reactant_item(world_->choice(rule_.reactants().at(0)));
    const ReactionRule::product_container_type& products(rule_.products());

    switch (products.size())
    {
        case 0:
            {
                world_->remove_voxel(reactant_item.voxel);
                ReactionInfo rinfo(world_->t());
                rinfo.add_reactant(reactant_item);
                push_reaction(std::make_pair(rule_, rinfo));
            }
            break;
        case 1:
            push_reaction(std::make_pair(rule_,
                                         apply_a2b(world_, reactant_item, *(products.begin()))));
            break;
        case 2:
            {
                ReactionInfo rinfo(apply_a2bc(world_, reactant_item,
                            *(products.begin()), (*(++products.begin()))));
                if (rinfo.has_occurred())
                    push_reaction(std::make_pair(rule_, rinfo));
            }
            break;
    }
    time_ += draw_dt();
}

Real FirstOrderReactionEvent::draw_dt()
{
    const Species& reactant(*(rule_.reactants().begin()));
    const Integer num_r(world_->num_voxels_exact(reactant));
    const Real k(rule_.k());
    const Real p = k * num_r;
    Real dt(inf);
    if (p != 0.)
    {
        const Real rnd(world_->rng()->uniform(0.,1.));
        dt = - log(1 - rnd) / p;
    }
    return dt;
}

} // spatiocyte

} // ecell4
