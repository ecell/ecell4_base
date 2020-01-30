#include "SpatiocyteEvent.hpp"
#include "utils.hpp"

namespace ecell4
{

namespace spatiocyte
{

StepEvent::StepEvent(boost::shared_ptr<Model> model,
                     boost::shared_ptr<SpatiocyteWorld> world,
                     const Species &species, const Real &t, const Real alpha)
    : SpatiocyteEvent(t), model_(model), world_(world), alpha_(alpha)
{
    if (const auto space_and_molecule_pool =
            world_->find_space_and_molecule_pool(species))
    {
        space_ = space_and_molecule_pool->first;
        mpool_ = space_and_molecule_pool->second;
    }
    else
    {
        throw "MoleculePool is not found";
    }

    time_ = t;
}

StepEvent3D::StepEvent3D(boost::shared_ptr<Model> model,
                         boost::shared_ptr<SpatiocyteWorld> world,
                         const Species &species, const Real &t,
                         const Real alpha)
    : StepEvent(model, world, species, t, alpha)
{
    const MoleculeInfo minfo(world_->get_molecule_info(species));
    const Real D(minfo.D);
    const Real R(world_->voxel_radius());

    if (D <= 0)
        dt_ = inf;
    else
        dt_ = 2 * R * R / 3 / D * alpha_;

    time_ = t + dt_;
}

StepEvent2D::StepEvent2D(boost::shared_ptr<Model> model,
                         boost::shared_ptr<SpatiocyteWorld> world,
                         const Species &species, const Real &t,
                         const Real alpha)
    : StepEvent(model, world, species, t, alpha)
{
    const MoleculeInfo minfo(world_->get_molecule_info(species));
    const Real D(minfo.D);
    const Real R(world_->voxel_radius());

    if (D <= 0)
        dt_ = inf;
    else
        dt_ = R * R / D * alpha_;

    time_ = t + dt_;
}

void StepEvent::attempt_reaction_(
    const SpatiocyteWorld::coordinate_id_pair_type &info, const Voxel &dst,
    const Real &alpha)
{
    const Voxel voxel(space_, info.coordinate);
    boost::shared_ptr<const VoxelPool> from_mt(voxel.get_voxel_pool());
    boost::shared_ptr<const VoxelPool> to_mt(dst.get_voxel_pool());

    const Species &speciesA(from_mt->species());
    const Species &speciesB(to_mt->species());

    const std::vector<ReactionRule> rules(
        model_->query_reaction_rules(speciesA, speciesB));

    if (rules.empty())
    {
        return;
    }

    const Real from_D(world_->get_molecule_info(speciesA).D);
    const Real to_D(world_->get_molecule_info(speciesB).D);
    const Real factor(
        calculate_dimensional_factor(from_mt, from_D, to_mt, to_D, world_));
    const Real rnd(world_->rng()->uniform(0, 1));
    Real accp(0.0);

    for (const auto &rule : rules)
    {
        const Real k(rule.k());
        const Real P(k * factor * alpha);
        accp += P;
        if (accp > 1 && k != std::numeric_limits<Real>::infinity())
        {
            std::cerr << "The total acceptance probability [" << accp
                      << "] exceeds 1 for '" << speciesA.serial() << "' and '"
                      << speciesB.serial() << "'." << std::endl;
        }
        if (accp >= rnd)
        {
            ReactionInfo rinfo(apply_second_order_reaction(
                world_, rule,
                ReactionInfo::Item(info.pid, from_mt->species(), voxel),
                ReactionInfo::Item(to_mt->get_particle_id(dst.coordinate),
                                   to_mt->species(), dst)));
            if (rinfo.has_occurred())
            {
                reaction_type reaction(std::make_pair(rule, rinfo));
                push_reaction(reaction);
            }
            return;
        }
    }
}

} // namespace spatiocyte

} // namespace ecell4
