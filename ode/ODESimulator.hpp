#ifndef __ECELL4_ODE_ODE_SIMULATOR_HPP
#define __ECELL4_ODE_ODE_SIMULATOR_HPP

#include <cstring>
#include <vector>
#include <boost/shared_ptr.hpp>

#include <ecell4/core/exceptions.hpp>
#include <ecell4/core/types.hpp>
#include <ecell4/core/get_mapper_mf.hpp>
#include <ecell4/core/NetworkModel.hpp>
#include <ecell4/core/Simulator.hpp>

#include "ODEWorld.hpp"

#include <boost/scoped_ptr.hpp>
#include <boost/scoped_array.hpp>

namespace ecell4
{

namespace ode
{

class ODESystem
{
public:

    typedef std::vector<double> state_type;

protected:

    typedef utils::get_mapper_mf<
        Species, state_type::size_type>::type species_map_type;

public:

    ODESystem(boost::shared_ptr<NetworkModel> model, const Real& volume)
        : model_(model), volume_(volume)
    {
        initialize();
    }

    void initialize()
    {
        const std::vector<Species> species(model_->list_species());
        state_type::size_type i(0);
        for (std::vector<Species>::const_iterator
                 it(species.begin()); it != species.end(); ++it)
        {
            index_map_[*it] = i;
            ++i;
        }
    }

    void operator()(const state_type& x, state_type& dxdt, const double& t)
    {
        for (state_type::iterator i(dxdt.begin()); i != dxdt.end(); ++i)
        {
            *i = 0.0;
        }

        const NetworkModel::reaction_rule_container_type&
            reaction_rules(model_->reaction_rules());
        for (NetworkModel::reaction_rule_container_type::const_iterator
                 i(reaction_rules.begin()); i != reaction_rules.end(); ++i)
        {
            double flux((*i).k() * volume_);

            const ReactionRule::reactant_container_type&
                reactants((*i).reactants());
            const ReactionRule::product_container_type&
                products((*i).products());
            for (ReactionRule::reactant_container_type::iterator
                     j(reactants.begin()); j != reactants.end(); ++j)
            {
                flux *= x[index_map_[*j]] / volume_;
            }

            for (ReactionRule::reactant_container_type::iterator
                     j(reactants.begin()); j != reactants.end(); ++j)
            {
                dxdt[index_map_[*j]] -= flux;
            }

            for (ReactionRule::product_container_type::iterator
                     j(products.begin()); j != products.end(); ++j)
            {
                dxdt[index_map_[*j]] += flux;
            }
        }
    }

protected:

    boost::shared_ptr<NetworkModel> model_;
    Real volume_;

    species_map_type index_map_;
};

struct StateAndTimeBackInserter
{
    typedef std::vector<ODESystem::state_type> state_container_type;
    typedef std::vector<double> time_container_type;

    state_container_type& m_states;
    time_container_type& m_times;

    StateAndTimeBackInserter(
        state_container_type& states, time_container_type& times)
        : m_states(states), m_times(times)
    {
        ;
    }

    void operator()(const ODESystem::state_type&x, double t)
    {
        m_states.push_back(x);
        m_times.push_back(t);
    }
};

class ODESimulator
    : public Simulator
{
public:

    ODESimulator(
        boost::shared_ptr<NetworkModel> model,
        boost::shared_ptr<ODEWorld> world)
        : model_(model), world_(world), dt_(0.0), num_steps_(0), is_dirty_(true)
    {
        ;
    }

    void initialize()
    {
        if (!is_dirty_)
        {
            return;
        }

        const std::vector<Species> species((*model_).list_species());
        for (std::vector<Species>::const_iterator
                 i(species.begin()); i != species.end(); ++i)
        {
            if (!(*world_).has_species(*i))
            {
                (*world_).reserve_species(*i);
            }
        }

        is_dirty_ = false;
    }

    // SimulatorTraits

    Real t(void) const
    {
        return (*world_).t();
    }

    Integer num_steps(void) const
    {
        return num_steps_;
    }

    Real dt() const
    {
        return dt_;
    }

    void step(void)
    {
        step(next_time());
    }

    bool step(const Real& upto);

    // Optional members

    void set_t(const Real& t)
    {
        (*world_).set_t(t);
    }

    void set_dt(const Real& dt)
    {
        if (dt <= 0)
        {
            throw std::invalid_argument("The step size must be positive.");
        }
        dt_ = dt;
    }

protected:

    boost::shared_ptr<NetworkModel> model_;
    boost::shared_ptr<ODEWorld> world_;
    Real dt_;
    Integer num_steps_;
    bool is_dirty_;
};

} // ode

} // ecell4

#endif /* __ECELL4_ODE_ODE_SIMULATOR_HPP */
