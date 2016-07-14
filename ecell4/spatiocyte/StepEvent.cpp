#include "SpatiocyteEvent.hpp"
#include "utils.hpp"

namespace ecell4
{

namespace spatiocyte
{

StepEvent::StepEvent(boost::shared_ptr<Model> model, boost::shared_ptr<SpatiocyteWorld> world,
        const Species& species, const Real& t, const Real alpha)
    : SpatiocyteEvent(t), model_(model), world_(world), species_(species), alpha_(alpha)
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

void StepEvent::fire_()
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
            attempt_reaction_(info, neighbor, alpha);
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
            else if (attempt_reaction_(info, neighbor, alpha).first != NO_REACTION)
                break;
        }
        ++idx;
    }
}

std::pair<StepEvent::attempt_reaction_result_type, StepEvent::reaction_type>
StepEvent::attempt_reaction_(
    const SpatiocyteWorld::coordinate_id_pair_type& info,
    const SpatiocyteWorld::coordinate_type to_coord,
    const Real& alpha)
{
    const VoxelPool* from_mt(
        world_->find_voxel_pool(info.coordinate));
    const VoxelPool* to_mt(
        world_->find_voxel_pool(to_coord));

    if (to_mt->is_vacant())
    {
        return std::make_pair(NO_REACTION, reaction_type());
    }

    const Species& speciesA(from_mt->species());
    const Species& speciesB(to_mt->species());

    const std::vector<ReactionRule> rules(
        model_->query_reaction_rules(speciesA, speciesB));

    if (rules.empty())
    {
        return std::make_pair(NO_REACTION, reaction_type());
    }

    const Real factor(calculate_dimensional_factor(from_mt, to_mt,
                boost::const_pointer_cast<const SpatiocyteWorld>(world_)));

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
                        world_->make_pid_voxel_pair(from_mt, info),
                        world_->make_pid_voxel_pair(to_mt, to_coord)));
            if (rinfo.has_occurred())
            {
                reaction_type reaction(std::make_pair(*itr, rinfo));
                push_reaction(reaction);
                return std::make_pair(REACTION_SUCCEEDED, reaction);
            }
            return std::make_pair(REACTION_FAILED, std::make_pair(*itr, rinfo));
        }
    }
    return std::make_pair(REACTION_FAILED, reaction_type());
}

} // spatiocyte

} // ecell4
