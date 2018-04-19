#include "SpatiocyteEvent.hpp"
#include "utils.hpp"

namespace ecell4
{

namespace spatiocyte
{

StepEvent::StepEvent(boost::shared_ptr<Model> model, boost::shared_ptr<SpatiocyteWorld> world,
        const Species& species, const Real& t, const Real alpha)
    : SpatiocyteEvent(t),
      model_(model),
      world_(world),
      mpool_(world_->find_molecule_pool(species)),
      alpha_(alpha)
{
    time_ = t;
}

StepEvent3D::StepEvent3D(boost::shared_ptr<Model> model,
                         boost::shared_ptr<SpatiocyteWorld> world,
                         const Species& species,
                         const Real& t,
                         const Real alpha)
    : StepEvent(model, world, species, t, alpha)
{
    const SpatiocyteWorld::molecule_info_type minfo(world_->get_molecule_info(species));
    const Real D(minfo.D);
    const Real R(world_->voxel_radius());

    if (D <= 0)
        dt_ = inf;
    else
        dt_ = 2 * R * R / 3 / D * alpha_;

    time_ = t + dt_;
}

void StepEvent3D::walk(const Real& alpha)
{
    if (alpha < 0 || alpha > 1)
    {
        return; // INVALID ALPHA VALUE
    }

    const boost::shared_ptr<RandomNumberGenerator>& rng(world_->rng());
    MoleculePool::container_type voxels;
    copy(mpool_->begin(), mpool_->end(), back_inserter(voxels));

    std::size_t idx(0);
    for (MoleculePool::container_type::iterator itr(voxels.begin());
         itr != voxels.end(); ++itr)
    {
        const SpatiocyteWorld::coordinate_id_pair_type& info(*itr);
        const Integer rnd(rng->uniform_int(0, world_->num_neighbors(info.coordinate)-1));

        if (world_->get_voxel_pool_at(info.coordinate) != mpool_)
        {
            // should skip if a voxel is not the target species.
            // when reaction has occured before, a voxel can be changed.
            continue;
        }

        const SpatiocyteWorld::coordinate_type
            neighbor(world_->get_neighbor(info.coordinate, rnd));

        if (world_->can_move(info.coordinate, neighbor))
        {
            if (rng->uniform(0,1) <= alpha)
                world_->move(info.coordinate, neighbor, /*candidate=*/idx);
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
                         const Species& species,
                         const Real& t,
                         const Real alpha)
    : StepEvent(model, world, species, t, alpha)
{
    const SpatiocyteWorld::molecule_info_type minfo(world_->get_molecule_info(species));
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

void StepEvent2D::walk(const Real& alpha)
{
    if (alpha < 0 || alpha > 1)
    {
        return; // INVALID ALPHA VALUE
    }

    const boost::shared_ptr<RandomNumberGenerator>& rng(world_->rng());
    MoleculePool::container_type voxels;
    copy(mpool_->begin(), mpool_->end(), back_inserter(voxels));

    std::size_t idx(0);
    for (MoleculePool::container_type::iterator itr(voxels.begin());
         itr != voxels.end(); ++itr)
    {
        const SpatiocyteWorld::coordinate_id_pair_type& info(*itr);
        if (world_->get_voxel_pool_at(info.coordinate) != mpool_)
        {
            // should skip if a voxel is not the target species.
            // when reaction has occured before, a voxel can be changed.
            continue;
        }

        const std::size_t num_neighbors(world_->num_neighbors(info.coordinate));

        ecell4::shuffle(*(rng.get()), nids_);
        for (std::vector<unsigned int>::const_iterator itr(nids_.begin());
             itr != nids_.end(); ++itr)
        {
            if (*itr >= num_neighbors)
                continue;

            const SpatiocyteWorld::coordinate_type neighbor(
                    world_->get_neighbor(info.coordinate, *itr));
            boost::shared_ptr<const VoxelPool> target(world_->get_voxel_pool_at(neighbor));

            if (target->get_dimension() > mpool_->get_dimension())
                continue;

            if (world_->can_move(info.coordinate, neighbor))
            {
                if (rng->uniform(0,1) <= alpha)
                    world_->move(info.coordinate, neighbor, /*candidate=*/idx);
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
    const SpatiocyteWorld::coordinate_id_pair_type& info,
    const SpatiocyteWorld::coordinate_type to_coord,
    const Real& alpha)
{
    boost::shared_ptr<const VoxelPool> from_mt(world_->get_voxel_pool_at(info.coordinate));
    boost::shared_ptr<const VoxelPool> to_mt(world_->get_voxel_pool_at(to_coord));

    if (to_mt->is_vacant())
    {
        return;
    }

    const Species& speciesA(from_mt->species());
    const Species& speciesB(to_mt->species());

    const std::vector<ReactionRule> rules(model_->query_reaction_rules(speciesA, speciesB));

    if (rules.empty())
    {
        return;
    }

    const Real factor(calculate_dimensional_factor(from_mt, to_mt, world_));
    const Real rnd(world_->rng()->uniform(0,1));
    Real accp(0.0);

    for (std::vector<ReactionRule>::const_iterator itr(rules.begin()); itr != rules.end(); ++itr)
    {
        const Real k((*itr).k());
        const Real P(k * factor * alpha);
        accp += P;
        if (accp > 1)
        {
            std::cerr << "The total acceptance probability [" << accp
                << "] exceeds 1 for '" << speciesA.serial()
                << "' and '" << speciesB.serial() << "'." << std::endl;
        }
        if (accp >= rnd)
        {
            ReactionInfo rinfo(apply_second_order_reaction(
                        world_, *itr,
                        ReactionInfo::Item(info.pid, from_mt->species(), info.coordinate),
                        ReactionInfo::Item(to_mt->get_particle_id(to_coord),
                                           to_mt->species(), to_coord)));
            if (rinfo.has_occurred())
            {
                reaction_type reaction(std::make_pair(*itr, rinfo));
                push_reaction(reaction);
            }
            return;
        }
    }
}

} // spatiocyte

} // ecell4
