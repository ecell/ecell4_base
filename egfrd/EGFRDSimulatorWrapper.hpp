#ifndef __EGFRD_SIMULATOR_HPP
#define __EGFRD_SIMULATOR_HPP

#include <stdexcept>
#include <boost/shared_ptr.hpp>

// epdp
#include "SpeciesTypeID.hpp"
#include "ReactionRule.hpp"
#include "EGFRDSimulator.hpp"
// epdp

#include <ecell4/core/NetworkModel.hpp>
#include <ecell4/core/Simulator.hpp>

#include "EGFRDWorld.hpp"


namespace ecell4
{

namespace egfrd
{

class EGFRDSimulatorWrapper
    : public Simulator
{
public:

    typedef EGFRDWorld::simulator_type simulator_type;

protected:

    typedef simulator_type::world_type world_type;
    typedef simulator_type::traits_type::network_rules_type network_rules_type;

public:

    EGFRDSimulatorWrapper(
        boost::shared_ptr<NetworkModel> model,
        boost::shared_ptr<EGFRDWorld> world)
        : model_(model), world_(world), rng_(),
          num_steps_(0), dt_(0)
    {
        NetworkModel::species_container_type const&
            species((*model_).species());
        for (NetworkModel::species_container_type::const_iterator
                 i(species.begin()); i != species.end(); ++i)
        {
            (*world_).add_species(*i);
        }

        NetworkModel::reaction_rule_container_type const&
            reaction_rules((*model_).reaction_rules());
        for (NetworkModel::reaction_rule_container_type::const_iterator
                 i(reaction_rules.begin()); i != reaction_rules.end(); ++i)
        {
            (*world_).add_reaction_rule(*i);
        }

        int dissociation_retry_moves(3);
        sim_ = boost::shared_ptr<simulator_type>(
            new simulator_type(
                world->world(), boost::shared_ptr<network_rules_type>(
                    new network_rules_type(world->model().network_rules())),
                rng_, dissociation_retry_moves));
    }

    Real t() const
    {
        return (*sim_).t();
    }

    void set_t(Real const& t)
    {
        throw NotImplemented("Not implemented yet.");
    }

    Real dt() const
    {
        return (*sim_).dt();
    }

    void set_dt(Real const& dt)
    {
        throw NotSupported("Not supported.");
    }

    Integer num_steps() const
    {
        return (*sim_).num_steps();
    }

    void step();
    bool step(Real const& upto);

protected:

    boost::shared_ptr<NetworkModel> model_;
    boost::shared_ptr<EGFRDWorld> world_;
    world_type::traits_type::rng_type rng_;
    boost::shared_ptr<simulator_type> sim_;

    /**
     * the protected internal state of EGFRDSimulator.
     * they are needed to be saved/loaded with Visitor pattern.
     */
    Real dt_;
    Integer num_steps_;
};

} // egfrd

} // ecell4

#endif /* __EGFRD_SIMULATOR_HPP */
