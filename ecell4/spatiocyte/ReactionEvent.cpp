#include "SpatiocyteEvent.hpp"

namespace ecell4
{

namespace spatiocyte
{

/// ZerothOrderReactionEvent

ZerothOrderReactionEvent::ZerothOrderReactionEvent(
    boost::shared_ptr<SpatiocyteWorld> world, const ReactionRule &rule,
    const Real &t)
    : SpatiocyteEvent(t), world_(world), rule_(rule)
{
    time_ = t + draw_dt();
}

void ZerothOrderReactionEvent::fire_()
{
    ReactionInfo rinfo(world_->t());

    for (const auto &sp : rule_.products())
    {
        const MoleculeInfo info(world_->get_molecule_info(sp));

        if (const auto space_and_location =
                world_->find_space_and_voxel_pool(Species(info.loc)))
        {
            const auto space = space_and_location->first;
            const auto location = space_and_location->second;

            if (location->size() == 0)
            {
                time_ += draw_dt();
                return;
            }

            while (true)
            {
                const Voxel voxel(
                    space, world_->rng()->uniform_int(0, space->size() - 1));

                if (voxel.get_voxel_pool() != location)
                {
                    continue;
                }

                if (boost::optional<ParticleID> pid =
                        world_->new_particle(sp, voxel))
                {
                    rinfo.add_product(ReactionInfo::Item(*pid, sp, voxel));
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
    const Real p = k * world_->volume();
    Real dt(inf);
    if (p != 0.)
    {
        const Real rnd(world_->rng()->uniform(0., 1.));
        dt = -log(1 - rnd) / p;
    }
    return dt;
}

/// FirstOrderReactionEvent

FirstOrderReactionEvent::FirstOrderReactionEvent(
    boost::shared_ptr<SpatiocyteWorld> world, const ReactionRule &rule,
    const Real &t)
    : SpatiocyteEvent(t), world_(world), rng_(world->rng()), rule_(rule)
{
    // assert(rule_.reactants().size() == 1);
    time_ = t + draw_dt();
}

void FirstOrderReactionEvent::fire_()
{
    const ReactionInfo::Item reactant_item(choice());
    const ReactionRule::product_container_type &products(rule_.products());

    switch (products.size())
    {
    case 0: {
        reactant_item.voxel.clear();
        ReactionInfo rinfo(world_->t());
        rinfo.add_reactant(reactant_item);
        push_reaction(std::make_pair(rule_, rinfo));
    }
    break;
    case 1:
        push_reaction(std::make_pair(
            rule_, apply_a2b(world_, reactant_item, *(products.begin()))));
        break;
    case 2: {
        ReactionInfo rinfo(apply_a2bc(world_, reactant_item,
                                      *(products.begin()),
                                      (*(++products.begin()))));
        if (rinfo.has_occurred())
            push_reaction(std::make_pair(rule_, rinfo));
    }
    break;
    }
    time_ += draw_dt();
}

Real FirstOrderReactionEvent::draw_dt()
{
    const Species &reactant(*(rule_.reactants().begin()));
    const Integer num_r(world_->num_voxels_exact(reactant));
    const Real k(rule_.k());
    if (num_r > 0)
    {
        const Real p = k * num_r;
        Real dt(inf);
        if (p != 0.)
        {
            const Real rnd(world_->rng()->uniform(0., 1.));
            dt = -log(1 - rnd) / p;
        }
        return dt;
    }
    else
    {
        return inf;
    }
}

} // namespace spatiocyte

} // namespace ecell4
