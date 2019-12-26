#include "SpatiocyteEvent.hpp"
#include "utils.hpp"

namespace ecell4
{

namespace spatiocyte
{

StepEvent::StepEvent(boost::shared_ptr<Model> model,
                     boost::shared_ptr<SpatiocyteWorld> world,
                     const Species &species, const Real &t, const Real alpha)
    : SpatiocyteEvent(t), model_(model), world_(world),
      mpool_(world_->find_molecule_pool(species)), alpha_(alpha)
{
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

void StepEvent3D::walk(const Real &alpha)
{
    if (alpha < 0 || alpha > 1)
    {
        return; // INVALID ALPHA VALUE
    }

    const boost::shared_ptr<RandomNumberGenerator> &rng(world_->rng());
    MoleculePool::container_type voxels;
    copy(mpool_->begin(), mpool_->end(), back_inserter(voxels));

    std::size_t idx(0);
    for (const auto &info : voxels)
    {
        const Voxel voxel(world_->coordinate2voxel(info.coordinate));

        if (voxel.get_voxel_pool() != mpool_)
        {
            // should skip if a voxel is not the target species.
            // when reaction has occured before, a voxel can be changed.
            continue;
        }

        const Voxel neighbor(voxel.get_neighbor_randomly(rng));

        if (world_->can_move(voxel, neighbor))
        {
            if (rng->uniform(0, 1) <= alpha)
                world_->move(voxel, neighbor, /*candidate=*/idx);
        }
        else
        {
            attempt_reaction_(info, neighbor, alpha);
        }

        ++idx;
    }
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

    nids_.clear();
    for (unsigned int i(0); i < 12; ++i)
        nids_.push_back(i);
}

void StepEvent2D::walk(const Real &alpha)
{
    if (alpha < 0 || alpha > 1)
    {
        return; // INVALID ALPHA VALUE
    }

    const boost::shared_ptr<RandomNumberGenerator> &rng(world_->rng());
    MoleculePool::container_type voxels;
    copy(mpool_->begin(), mpool_->end(), back_inserter(voxels));

    std::size_t idx(0);
    for (const auto &info : voxels)
    {
        // TODO: Calling coordinate2voxel is invalid
        const Voxel voxel(world_->coordinate2voxel(info.coordinate));

        if (voxel.get_voxel_pool() != mpool_)
        {
            // should skip if a voxel is not the target species.
            // when reaction has occured before, a voxel can be changed.
            continue;
        }

        const std::size_t num_neighbors(voxel.num_neighbors());

        ecell4::shuffle(*(rng.get()), nids_);
        for (const auto &index : nids_)
        {
            if (index >= num_neighbors)
                continue;

            const Voxel neighbor(voxel.get_neighbor(index));
            boost::shared_ptr<const VoxelPool> target(
                neighbor.get_voxel_pool());

            if (world_->get_dimension(target->species()) > Shape::TWO)
                continue;

            if (world_->can_move(voxel, neighbor))
            {
                if (rng->uniform(0, 1) <= alpha)
                    world_->move(voxel, neighbor, /*candidate=*/idx);
            }
            else
            {
                attempt_reaction_(info, neighbor, alpha);
            }
            break;
        }
        ++idx;
    }
}

void StepEvent::attempt_reaction_(
    const SpatiocyteWorld::coordinate_id_pair_type &info, const Voxel &dst,
    const Real &alpha)
{
    // TODO: Calling coordiante2voxel is invalid
    const Voxel voxel(world_->coordinate2voxel(info.coordinate));
    boost::shared_ptr<const VoxelPool> from_mt(voxel.get_voxel_pool());
    boost::shared_ptr<const VoxelPool> to_mt(dst.get_voxel_pool());

    if (to_mt->is_vacant())
    {
        return;
    }

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
