#ifndef __ECELL4_LATTICE_LATTICE_SIMULATOR_HPP
#define __ECELL4_LATTICE_LATTICE_SIMULATOR_HPP

#include <numeric>
#include <boost/shared_ptr.hpp>
#include <boost/lexical_cast.hpp>

#include <ecell4/core/Model.hpp>
#include <ecell4/core/ReactionRule.hpp>
#include <ecell4/core/Reaction.hpp>
#include <ecell4/core/MolecularTypeBase.hpp>
#include <ecell4/core/SimulatorBase.hpp>
#include <ecell4/core/RandomNumberGenerator.hpp>
#include <ecell4/core/EventScheduler.hpp>

#include "LatticeWorld.hpp"

namespace ecell4
{

namespace lattice
{

class LatticeSimulator
    : public SimulatorBase<Model, LatticeWorld>
{
public:

    typedef SimulatorBase<Model, LatticeWorld> base_type;
    typedef Reaction<Voxel> reaction_type;

protected:

    struct StepEvent : EventScheduler::Event
    {
        StepEvent(LatticeSimulator* sim, const Species& species, const Real& t)
            : EventScheduler::Event(t), sim_(sim), species_(species), alpha_(1.0)
        {
            const LatticeWorld::molecule_info_type
                minfo(sim_->world_->get_molecule_info(species));
            const Real R(minfo.radius);
            const Real D(minfo.D);
            // const Real R(sim_->world_->voxel_radius());
            // Real D = boost::lexical_cast<Real>(species.get_attribute("D"));
            if (D <= 0)
            {
                dt_ = inf;
            } else {
                dt_ = 2 * R * R / 3 / D * alpha_;
            }

            time_ = t + dt_;
        }

        virtual ~StepEvent()
        {
        }

        virtual void fire()
        {
            sim_->walk(species_, alpha_);
            time_ += dt_;
        }

        Species const& species() const
        {
            return species_;
        }

        Real const& alpha() const
        {
            return alpha_;
        }

    protected:

        LatticeSimulator* sim_;
        Species species_;
        MolecularTypeBase* mt_;
        Real alpha_;
    };

    struct ZerothOrderReactionEvent : EventScheduler::Event
    {
        ZerothOrderReactionEvent(
            LatticeSimulator* sim, const ReactionRule& rule, const Real& t)
            : EventScheduler::Event(t), sim_(sim), rule_(rule)
        {
            time_ = t + draw_dt();
        }

        virtual ~ZerothOrderReactionEvent()
        {
        }

        virtual void fire()
        {
            sim_->apply_zeroth_order_reaction_(rule_);
            time_ += draw_dt();
        }

        virtual void interrupt(Real const& t)
        {
            time_ = t + draw_dt();
        }

        Real draw_dt()
        {
            const Real k(rule_.k());
            const Real p = k;
            Real dt(inf);
            if (p != 0.)
            {
                const Real rnd(sim_->world_->rng()->uniform(0.,1.));
                dt = - log(1 - rnd) / p;
            }
            return dt;
        }

    protected:

        LatticeSimulator* sim_;
        ReactionRule rule_;
    };

    struct FirstOrderReactionEvent : EventScheduler::Event
    {
        FirstOrderReactionEvent(
            LatticeSimulator* sim, const ReactionRule& rule, const Real& t)
            : EventScheduler::Event(t), sim_(sim), rule_(rule)
        {
            time_ = t + draw_dt();
        }

        virtual ~FirstOrderReactionEvent()
        {
        }

        virtual void fire()
        {
            const Species reactant(*(rule_.reactants().begin()));
            sim_->apply_first_order_reaction_(rule_, sim_->world_->choice(reactant));
            time_ += draw_dt();
        }

        virtual void interrupt(Real const& t)
        {
            time_ = t + draw_dt();
        }

        Real draw_dt()
        {
            const Species reactant(*(rule_.reactants().begin()));
            const Integer num_r(sim_->world_->num_voxels_exact(reactant));
            const Real k(rule_.k());
            const Real p = k * num_r;
            Real dt(inf);
            if (p != 0.)
            {
                const Real rnd(sim_->world_->rng()->uniform(0.,1.));
                dt = - log(1 - rnd) / p;
            }
            return dt;
        }

    protected:

        LatticeSimulator* sim_;
        ReactionRule rule_;
    };

public:

    LatticeSimulator(
            boost::shared_ptr<Model> model,
            boost::shared_ptr<LatticeWorld> world)
        : base_type(model, world)
    {
        initialize();
    }

    LatticeSimulator(boost::shared_ptr<LatticeWorld> world)
        : base_type(world)
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
    // void run(const Real& duration);
    void walk(const Species& species);
    void walk(const Species& species, const Real& alpha);

    std::vector<ReactionRule> last_reactions() const
    {
        return reactions_;
    }

protected:

    boost::shared_ptr<EventScheduler::Event> create_step_event(
        const Species& species, const Real& t);
    boost::shared_ptr<EventScheduler::Event> create_zeroth_order_reaction_event(
        const ReactionRule& reaction_rule, const Real& t);
    boost::shared_ptr<EventScheduler::Event> create_first_order_reaction_event(
        const ReactionRule& reaction_rule, const Real& t);
    std::pair<bool, reaction_type> attempt_reaction_(
        const LatticeWorld::particle_info_type info,
        LatticeWorld::coordinate_type to_coord);
    std::pair<bool, reaction_type> apply_second_order_reaction_(
        const ReactionRule& reaction_rule,
        const reaction_type::particle_type& p0,
        const reaction_type::particle_type& p1);
    std::pair<bool, reaction_type> apply_first_order_reaction_(
        const ReactionRule& reaction_rule,
        const reaction_type::particle_type& p);
    void apply_ab2c(
        const LatticeWorld::private_coordinate_type coord,
        const Species& product_species,
        reaction_type& reaction);
    void apply_ab2cd(
        const LatticeWorld::private_coordinate_type from_coord,
        const LatticeWorld::private_coordinate_type to_coord,
        const Species& product_species0,
        const Species& product_species1,
        reaction_type& reaction);
    void apply_a2b(
        const LatticeWorld::private_coordinate_type coord,
        const Species& product_species,
        reaction_type& reaction);
    bool apply_a2bc(
        const LatticeWorld::private_coordinate_type coord,
        const Species& product_species0,
        const Species& product_species1,
        reaction_type& reaction);
    std::pair<bool, reaction_type> apply_zeroth_order_reaction_(
        const ReactionRule& reaction_rule);

    void  register_product_species(const Species& product_species);

    void step_();
    void register_events(const Species& species);
    // void register_step_event(const Species& species);

    inline Voxel private_voxel2voxel(const Voxel& v) const
    {
        const LatticeWorld::coordinate_type
            coord(world_->private2coord(v.coordinate()));
        return Voxel(v.species(), coord, v.radius(), v.D(), v.loc());
    }

protected:

    EventScheduler scheduler_;
    std::vector<ReactionRule> reactions_;
    std::vector<Species> new_species_;

    Real dt_;
};

} // lattice

} // ecell4

#endif /* __ECELL4_LATTICE_LATTICE_SIMULATOR_HPP */
