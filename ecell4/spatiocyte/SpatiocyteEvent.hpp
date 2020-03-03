#ifndef ECELL4_SPATIOCYTE_EVENT_HPP
#define ECELL4_SPATIOCYTE_EVENT_HPP

#include "SpatiocyteReactions.hpp"
#include "SpatiocyteWorld.hpp"
#include "utils.hpp"
#include <ecell4/core/EventScheduler.hpp>
#include <ecell4/core/Model.hpp>
#include <ecell4/core/ReactionRule.hpp>

namespace ecell4
{

namespace spatiocyte
{

struct SpatiocyteEvent : public Event
{
public:
    typedef std::pair<ReactionRule, ReactionInfo> reaction_type;

    SpatiocyteEvent(Real const &time) : Event(time) {}
    virtual ~SpatiocyteEvent() {}

    const std::vector<reaction_type> &reactions() const { return reactions_; }

    virtual void fire()
    {
        reactions_.clear();
        fire_();
    }

protected:
    virtual void fire_() = 0;

    void push_reaction(const reaction_type &reaction)
    {
        reactions_.push_back(reaction);
    }

    std::vector<reaction_type> reactions_;
};

template <int Dimension>
const Real calc_dt(const Real R, const Real D);

template <int Dimension>
struct StepEvent : SpatiocyteEvent
{
    StepEvent(boost::shared_ptr<Model> model,
              boost::shared_ptr<SpatiocyteWorld> world, const Species &species,
              const Real &t, const Real alpha = 1.0)
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

        const MoleculeInfo minfo(world_->get_molecule_info(species));
        const Real D(minfo.D);
        const Real R(world_->voxel_radius());

        if (D <= 0)
            dt_ = std::numeric_limits<Real>::infinity();
        else
            dt_ = calc_dt<Dimension>(R, D) * alpha_;

        time_ = t + dt_;
    }

    Species const &species() const { return mpool_->species(); }

    Real const &alpha() const { return alpha_; }

    void fire_()
    {
        walk(alpha_);
        time_ += dt_;
    }

    void walk(const Real &alpha)
    {
        if (alpha < 0 || alpha > 1)
        {
            return; // INVALID ALPHA VALUE
        }

        MoleculePool::container_type voxels;
        copy(mpool_->begin(), mpool_->end(), back_inserter(voxels));

        std::size_t idx(0);
        for (const auto &info : voxels)
        {
            const Voxel voxel(space_, info.coordinate);

            if (voxel.get_voxel_pool() != mpool_)
            {
                // should skip if a voxel is not the target species.
                // when reaction has occured before, a voxel can be changed.
                continue;
            }

            const Voxel neighbor =
                world_->get_neighbor_randomly<Dimension>(voxel);

            if (world_->can_move(voxel, neighbor))
            {
                if (world_->rng()->uniform(0, 1) <= alpha)
                    world_->move(voxel, neighbor, /*candidate=*/idx);
            }
            else
            {
                attempt_reaction_(info, neighbor, alpha);
            }

            ++idx;
        }
    }

protected:
    void attempt_reaction_(const SpatiocyteWorld::coordinate_id_pair_type &info,
                           const Voxel &dst, const Real &alpha)
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
                          << "] exceeds 1 for '" << speciesA.serial()
                          << "' and '" << speciesB.serial() << "'."
                          << std::endl;
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

protected:
    boost::shared_ptr<Model> model_;
    boost::shared_ptr<SpatiocyteWorld> world_;
    boost::weak_ptr<VoxelSpaceBase> space_;
    boost::shared_ptr<MoleculePool> mpool_;

    const Real alpha_;
};

struct ZerothOrderReactionEvent : SpatiocyteEvent
{
    ZerothOrderReactionEvent(boost::shared_ptr<SpatiocyteWorld> world,
                             const ReactionRule &rule, const Real &t);

    virtual ~ZerothOrderReactionEvent() {}
    virtual void fire_();

    Real draw_dt();
    virtual void interrupt(Real const &t) { time_ = t + draw_dt(); }

protected:
    boost::shared_ptr<SpatiocyteWorld> world_;
    ReactionRule rule_;
};

struct FirstOrderReactionEvent : SpatiocyteEvent
{
    FirstOrderReactionEvent(boost::shared_ptr<SpatiocyteWorld> world,
                            const ReactionRule &rule, const Real &t);

    virtual ~FirstOrderReactionEvent() {}
    virtual void fire_();

    Real draw_dt();
    virtual void interrupt(Real const &t) { time_ = t + draw_dt(); }

protected:
    ReactionInfo::Item choice()
    {
        const Species &species(rule_.reactants().at(0));
        if (const auto space_and_molecule_pool =
                world_->find_space_and_molecule_pool(species))
        {
            const auto space = space_and_molecule_pool->first;
            const auto molecule_pool = space_and_molecule_pool->second;

            const auto i =
                rng_.lock()->uniform_int(0, molecule_pool->size() - 1);
            const auto &info = molecule_pool->at(i);

            return ReactionInfo::Item(info.pid, species,
                                      Voxel(space, info.coordinate));
        }
        throw "MoleculePool is not found";
    }

    boost::shared_ptr<SpatiocyteWorld> world_;
    boost::weak_ptr<RandomNumberGenerator> rng_;
    ReactionRule rule_;
};

} // namespace spatiocyte

} // namespace ecell4

#endif /* ECELL4_SPATIOCYTE_EVENT_HPP */
