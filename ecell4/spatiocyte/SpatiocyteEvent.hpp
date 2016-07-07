#ifndef __ECELL4_SPATIOCYTE_EVENT_HPP
#define __ECELL4_SPATIOCYTE_EVENT_HPP

#include <ecell4/core/ReactionRule.hpp>
#include <ecell4/core/EventScheduler.hpp>
#include <ecell4/core/VoxelPool.hpp>
#include "SpatiocyteWorld.hpp"

namespace ecell4
{

namespace spatiocyte {

class SpatiocyteSimulator;

class ReactionInfo
{
public:

    typedef std::pair<ParticleID, Voxel> particle_id_pair_type;
    typedef std::vector<particle_id_pair_type> container_type;

public:

    ReactionInfo() : t_(0), reactants_(), products_() {}

    ReactionInfo(const Real t) : t_(t), reactants_(), products_() {}

    ReactionInfo(
        const Real t,
        const container_type& reactants,
        const container_type& products)
        : t_(t), reactants_(reactants), products_(products) {}

    ReactionInfo(const ReactionInfo& another)
        : t_(another.t()), reactants_(another.reactants()), products_(another.products()) {}

    Real t() const
    {
        return t_;
    }

    const container_type& reactants() const
    {
        return reactants_;
    }

    void add_reactant(const particle_id_pair_type& pid_pair)
    {
        reactants_.push_back(pid_pair);
    }

    const container_type& products() const
    {
        return products_;
    }

    void add_product(const particle_id_pair_type& pid_pair)
    {
        products_.push_back(pid_pair);
    }

protected:

    Real t_;
    container_type reactants_, products_;
};

struct SpatiocyteEvent : public Event
{
public:
    typedef ReactionInfo reaction_info_type;
    typedef std::pair<ReactionRule, reaction_info_type> reaction_type;

    SpatiocyteEvent(Real const& time) : Event(time) {}
    virtual ~SpatiocyteEvent() {}

    std::vector<reaction_type> last_reactions() const
    {
        return last_reactions_;
    }

    void push_reaction(const reaction_type& reaction)
    {
        last_reactions_.push_back(reaction);
    }

protected:
    std::vector<reaction_type> last_reactions_;

};

struct StepEvent : SpatiocyteEvent
{
    StepEvent(SpatiocyteSimulator* sim, const Species& species, const Real& t,
        const Real alpha=1.0);

    virtual ~StepEvent() {}
    virtual void fire();

    Species const& species() const
    {
        return species_;
    }

    Real const& alpha() const
    {
        return alpha_;
    }

    void walk(const Real& alpha);

protected:

    void walk_in_space_(const MoleculePool* mtype, const Real& alpha);
    void walk_on_surface_(const MoleculePool* mtype, const Real& alpha);

    SpatiocyteSimulator* sim_;
    boost::shared_ptr<SpatiocyteWorld> world_;
    Species species_;
    VoxelPool* mt_;
    const Real alpha_;
    std::vector<unsigned int> nids_; // neighbor indexes
};

struct ZerothOrderReactionEvent : SpatiocyteEvent
{
    ZerothOrderReactionEvent(
        boost::shared_ptr<SpatiocyteWorld> world, const ReactionRule& rule, const Real& t);

    virtual ~ZerothOrderReactionEvent() {}
    virtual void fire();

    Real draw_dt();
    virtual void interrupt(Real const& t)
    {
        time_ = t + draw_dt();
    }

protected:

    boost::shared_ptr<SpatiocyteWorld> world_;
    ReactionRule rule_;
};

struct FirstOrderReactionEvent : SpatiocyteEvent
{
    FirstOrderReactionEvent(
        boost::shared_ptr<SpatiocyteWorld> world, const ReactionRule& rule, const Real& t);

    virtual ~FirstOrderReactionEvent() {}
    virtual void fire();

    Real draw_dt();
    virtual void interrupt(Real const& t)
    {
        time_ = t + draw_dt();
    }

protected:

    const std::string get_serial(const SpatiocyteWorld::coordinate_type coord) const
    {
        const VoxelPool* mtype(world_->find_voxel_pool(coord));
        return mtype->is_vacant() ? "" : mtype->species().serial();
    }

    const std::string get_location(const SpatiocyteWorld::coordinate_type coord) const
    {
        const VoxelPool* mtype(world_->find_voxel_pool(coord));
        if (mtype->is_vacant())
            return "";
        const VoxelPool* ltype(mtype->location());
        return ltype->is_vacant() ? "" : ltype->species().serial();
    }

    reaction_info_type apply_a2b(
        const reaction_info_type::particle_id_pair_type& p,
        const Species& product_species);

    std::pair<bool, reaction_info_type> apply_a2bc(
        const reaction_info_type::particle_id_pair_type& p,
        const Species& product_species0,
        const Species& product_species1);

    boost::shared_ptr<SpatiocyteWorld> world_;
    ReactionRule rule_;
};

} // spatiocyte

} // ecell4

#endif /* __ECELL4_SPATIOCYTE_EVENT_HPP */
