#ifndef __ECELL4_LATTICE_LATTICE_SIMULATOR_HPP
#define __ECELL4_LATTICE_LATTICE_SIMULATOR_HPP

#include <numeric>
#include <boost/shared_ptr.hpp>
#include <boost/lexical_cast.hpp>

#include <ecell4/core/Model.hpp>
#include <ecell4/core/ReactionRule.hpp>
// #include <ecell4/core/Reaction.hpp>
#include <ecell4/core/VoxelPool.hpp>
#include <ecell4/core/SimulatorBase.hpp>
#include <ecell4/core/RandomNumberGenerator.hpp>
#include <ecell4/core/EventScheduler.hpp>
#include <ecell4/core/get_mapper_mf.hpp>

#include "SpatiocyteWorld.hpp"
#include "SpatiocyteEvent.hpp"

namespace ecell4
{

namespace spatiocyte
{

class SpatiocyteSimulator
    : public SimulatorBase<Model, SpatiocyteWorld>
{
public:

    typedef SimulatorBase<Model, SpatiocyteWorld> base_type;

    typedef ReactionInfo reaction_info_type;
    typedef std::pair<ReactionRule, reaction_info_type> reaction_type;

    typedef EventSchedulerBase<SpatiocyteEvent> scheduler_type;

    typedef utils::get_mapper_mf<Species, Real>::type alpha_map_type;

public:

    SpatiocyteSimulator(
            boost::shared_ptr<Model> model,
            boost::shared_ptr<SpatiocyteWorld> world,
            const Real alpha = 1.0)
        : base_type(model, world), alpha_(alpha)
    {
        initialize();
    }

    SpatiocyteSimulator(
            boost::shared_ptr<SpatiocyteWorld> world,
            const Real alpha = 1.0)
        : base_type(world), alpha_(alpha)
    {
        initialize();
    }

    virtual Real dt() const
    {
        return dt_;
    }

    void initialize();
    void finalize();
    void step();
    bool step(const Real& upto);
    // void walk(const Species& species);
    void walk(const Species& species, const Real& alpha);
    Real calculate_alpha(const ReactionRule& rule) const;

    virtual bool check_reaction() const
    {
        return last_event_->last_reactions().size() > 0;
        // return last_reactions_.size() > 0;
    }

    std::vector<std::pair<ReactionRule, reaction_info_type> > last_reactions() const
    {
        return last_event_->last_reactions();
        // return last_reactions_;
    }

    void set_alpha(const Real alpha)
    {
        if (alpha < 0 || alpha > 1)
        {
            return;  // XXX: ValueError
        }

        alpha_ = alpha;
        initialize();
    }

    Real get_alpha() const
    {
        return alpha_;
    }

    // TODO: remove the below public functions
    boost::shared_ptr<SpatiocyteWorld> world()
    {
        return world_;
    }
    std::pair<bool, reaction_type> apply_zeroth_order_reaction_(
        const ReactionRule& reaction_rule);
    std::pair<bool, reaction_type> apply_first_order_reaction_(
        const ReactionRule& reaction_rule,
        const reaction_info_type::particle_id_pair_type& p);

protected:

    boost::shared_ptr<SpatiocyteEvent> create_step_event(
        const Species& species, const Real& t);
    boost::shared_ptr<SpatiocyteEvent> create_zeroth_order_reaction_event(
        const ReactionRule& reaction_rule, const Real& t);
    boost::shared_ptr<SpatiocyteEvent> create_first_order_reaction_event(
        const ReactionRule& reaction_rule, const Real& t);
    Real calculate_dimensional_factor(
        const VoxelPool* mt0, const VoxelPool* mt1) const;

    typedef enum
    {
        NO_REACTION = 0,
        REACTION_FAILED = 1,
        REACTION_SUCCEEDED = 2
    } attempt_reaction_result_type;

    std::pair<attempt_reaction_result_type, reaction_type> attempt_reaction_(
        const SpatiocyteWorld::coordinate_id_pair_type& info,
        SpatiocyteWorld::coordinate_type to_coord, const Real& alpha);


    std::pair<bool, reaction_type> apply_a2b(
        const ReactionRule& reaction_rule,
        const reaction_info_type::particle_id_pair_type& p,
        const Species& product_species);
    std::pair<bool, reaction_type> apply_a2bc(
        const ReactionRule& reaction_rule,
        const reaction_info_type::particle_id_pair_type& p,
        const Species& product_species0,
        const Species& product_species1);

    std::pair<bool, reaction_type> apply_second_order_reaction_(
        const ReactionRule& reaction_rule,
        const reaction_info_type::particle_id_pair_type& p0,
        const reaction_info_type::particle_id_pair_type& p1);
    std::pair<bool, reaction_type> apply_vanishment(
        const ReactionRule& reaction_rule,
        const reaction_info_type::particle_id_pair_type& p0,
        const reaction_info_type::particle_id_pair_type& p1);
    std::pair<bool, reaction_type> apply_ab2c(
        const ReactionRule& reaction_rule,
        const reaction_info_type::particle_id_pair_type& p0,
        const reaction_info_type::particle_id_pair_type& p1,
        const Species& product_species);
    std::pair<bool, reaction_type> apply_ab2cd(
        const ReactionRule& reaction_rule,
        const reaction_info_type::particle_id_pair_type& p0,
        const reaction_info_type::particle_id_pair_type& p1,
        const Species& product_species0,
        const Species& product_species1);
    std::pair<bool, reaction_type> apply_ab2cd_in_order(
        const ReactionRule& reaction_rule,
        const reaction_info_type::particle_id_pair_type& p0,
        const reaction_info_type::particle_id_pair_type& p1,
        const Species& product_species0,
        const Species& product_species1,
        const SpatiocyteWorld::coordinate_type coord0,
        const SpatiocyteWorld::coordinate_type coord1);

    void register_product_species(const Species& product_species);
    // void register_reactant_species(
    //     const SpatiocyteWorld::coordinate_id_pair_type pinfo, reaction_type& reaction) const;

    void step_();
    void register_events(const Species& species);
    // void register_step_event(const Species& species);
    void update_alpha_map();

    void walk_in_space_(const MoleculePool* mtype, const Real& alpha);
    void walk_on_surface_(const MoleculePool* mtype, const Real& alpha);

    const std::string get_serial(
        const SpatiocyteWorld::coordinate_type coord) const
    {
        const VoxelPool* mtype(world_->find_voxel_pool(coord));
        return mtype->is_vacant() ? "" : mtype->species().serial();
    }

    const std::string get_location(
        const SpatiocyteWorld::coordinate_type coord) const
    {
        const VoxelPool* mtype(world_->find_voxel_pool(coord));
        if (mtype->is_vacant())
        {
            return "";
        }
        const VoxelPool* ltype(mtype->location());
        return ltype->is_vacant() ? "" : ltype->species().serial();
    }

protected:

    scheduler_type scheduler_;
    std::vector<reaction_type> last_reactions_;
    boost::shared_ptr<SpatiocyteEvent> last_event_;
    std::vector<Species> new_species_;
    std::vector<unsigned int> nids_; // neighbor indexes
    alpha_map_type alpha_map_;
    //map<Species> alpha_map_;

    Real dt_;
    Real alpha_;
};

} // spatiocyte

} // ecell4

#endif /* __ECELL4_LATTICE_LATTICE_SIMULATOR_HPP */
