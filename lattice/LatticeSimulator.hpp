#ifndef __ECELL4_LATTICE_LATTICE_SIMULATOR_HPP
#define __ECELL4_LATTICE_LATTICE_SIMULATOR_HPP

#include <numeric>
#include <boost/shared_ptr.hpp>
#include <boost/lexical_cast.hpp>

#include <ecell4/core/NetworkModel.hpp>
#include <ecell4/core/ReactionRule.hpp>
#include <ecell4/core/Simulator.hpp>
#include <ecell4/core/RandomNumberGenerator.hpp>

#include "LatticeWorld.hpp"
#include "EventScheduler.hpp"

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
            const Real R(sim_->world_->normalized_voxel_radius());
            Real D = boost::lexical_cast<Real>(species.get_attribute("D"));
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
            const std::vector<Coord> coords(sim_->world_->list_coords(reactant));
            const Integer index(sim_->world_->rng()->uniform_int(0,coords.size() - 1));
            const Coord coord(coords.at(index));
            sim_->apply_reaction_(rule_, coord);
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
    void attempt_reaction_(Coord from_coord, Coord to_coord);
    void apply_reaction_(const ReactionRule& reaction_rule,
            const Coord& from_coord, const Coord& to_coord);
    void apply_reaction_(const ReactionRule& reaction_rule, const Coord& cood);
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
