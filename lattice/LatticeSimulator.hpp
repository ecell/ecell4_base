#ifndef __ECELL4_LATTICE_LATTICE_SIMULATOR_HPP
#define __ECELL4_LATTICE_LATTICE_SIMULATOR_HPP

#include <numeric>
#include <boost/shared_ptr.hpp>
#include <boost/lexical_cast.hpp>

#include <ecell4/core/NetworkModel.hpp>
#include <ecell4/core/ReactionRule.hpp>
#include <ecell4/core/Reaction.hpp>
#include <ecell4/core/MolecularTypeBase.hpp>
#include <ecell4/core/Simulator.hpp>
#include <ecell4/core/RandomNumberGenerator.hpp>
#include <ecell4/core/EventScheduler.hpp>

#include "LatticeWorld.hpp"

namespace ecell4
{

namespace lattice
{

class LatticeSimulator
    : public Simulator
{
protected:
    struct StepEvent : EventScheduler::Event
    {
        StepEvent(LatticeSimulator* sim, const Species& species, const Real& t)
            : EventScheduler::Event(t), sim_(sim), species_(species)
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
                dt_ = 2 * R * R / 3 / D;
            }

            time_ = t + dt_;
        }

        virtual ~StepEvent()
        {
        }

        virtual void fire()
        {
            sim_->walk(species_);
            time_ += dt_;
        }

    protected:
        LatticeSimulator* sim_;
        Species species_;
        MolecularTypeBase* mt_;
        Real dt_;
    };

    struct FirstOrderReactionEvent : EventScheduler::Event
    {
        FirstOrderReactionEvent(LatticeSimulator* sim, const ReactionRule& rule)
            : EventScheduler::Event(0.0), sim_(sim), rule_(rule)
        {
            time_ = sim_->t() + draw_dt();
        }

        virtual ~FirstOrderReactionEvent()
        {
        }

        virtual void fire()
        {
            const Species reactant(*(rule_.reactants().begin()));
            MolecularTypeBase* mt(sim_->world_->find_molecular_type(reactant));
            const Integer index(sim_->world_->rng()->uniform_int(0, mt->size() - 1));
            sim_->apply_reaction_(rule_, mt->at(index));
            time_ += draw_dt();
        }

        virtual void interrupt(Real const& t)
        {
            time_ = t + draw_dt();
        }

        Real draw_dt()
        {
            const Species reactant(*(rule_.reactants().begin()));
            const Integer num_r(sim_->world_->num_molecules(reactant));
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
            boost::shared_ptr<NetworkModel> model,
            boost::shared_ptr<LatticeWorld> world)
        : model_(model), world_(world), is_initialized_(false)
    {
        world_->bind_to(model);
    }

    virtual Real t() const
    {
        return world_->t();
    }

    virtual Real dt() const
    {
        // TODO
        return 0.;
    }

    Integer num_steps() const
    {
        // TODO
        return 0;
    }

    void initialize();
    void step();
    bool step(const Real& upto);
    void walk(const Species& species);

protected:
    boost::shared_ptr<EventScheduler::Event> create_step_event(
            const Species& species, const Real& t);
    boost::shared_ptr<EventScheduler::Event> create_first_order_reaction_event(
            const ReactionRule& reaction_rule);
    std::pair<bool, Reaction<Voxel> > attempt_reaction_(
            LatticeWorld::particle_info& info, LatticeWorld::coordinate_type to_coord);
    std::pair<bool, Reaction<Voxel> > apply_reaction_(
            const ReactionRule& reaction_rule, LatticeWorld::particle_info& from_info,
            const LatticeWorld::particle_info& to_info);
    std::pair<bool, Reaction<Voxel> > apply_reaction_(
            const ReactionRule& reaction_rule, LatticeWorld::particle_info& info);
    void step_();
    void register_step_event(const Species& species);

protected:

    boost::shared_ptr<NetworkModel> model_;
    boost::shared_ptr<LatticeWorld> world_;
    bool is_initialized_;

    EventScheduler scheduler_;
    std::vector<ReactionRule> reactions_;
    std::vector<Species> new_species_;

    Real dt_;

};

} // lattice

} // ecell4

#endif /* __ECELL4_LATTICE_LATTICE_SIMULATOR_HPP */
