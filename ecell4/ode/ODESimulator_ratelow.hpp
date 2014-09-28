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
#include <ecell4/core/ModelWrapper.hpp>

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
    typedef Model::reaction_rule_container_type reaction_rule_container_type;

protected:

    typedef utils::get_mapper_mf<
        Species, state_type::size_type>::type species_map_type;

public:

    ODESystem(const std::vector<Species>& species,
        const reaction_rule_container_type& reaction_rules, const Real& volume)
        : species_(species), reaction_rules_(reaction_rules), volume_(volume)
    {
        initialize();
    }

    void initialize()
    {
        state_type::size_type i(0);
        for (std::vector<Species>::const_iterator
                 it(species_.begin()); it != species_.end(); ++it)
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

        Ratelow::state_container_type reactant_states;
        reactant_states.reserve(5); // A little bigger than usual reactants' size of ReactionRules.
        for (reaction_rule_container_type::const_iterator
            i(reaction_rules_.begin()); i != reaction_rules_.end(); ++i)
        {
            const ReactionRule::reactant_container_type&
                reactants((*i).reactants());
            const ReactionRule::product_container_type&
                products((*i).products());

            //boost::scoped_array<Real> reactant_states(new Real[reactants.size()]);
            reactant_states.resize( reactants.size() );
            state_type::size_type cnt(0);
            for (ReactionRule::reactant_container_type::const_iterator
                     j(reactants.begin()); j != reactants.end(); ++j)
            {
                reactant_states[cnt] = x[index_map_[*j]];
            }
            boost::shared_ptr<Ratelow> ratelow = i->get_ratelow();
            double flux = (*ratelow)(reactant_states, volume_);
            for (ReactionRule::reactant_container_type::const_iterator
                     j(reactants.begin()); j != reactants.end(); ++j)
            {
                dxdt[index_map_[*j]] -= flux;
            }

            for (ReactionRule::product_container_type::const_iterator
                     j(products.begin()); j != products.end(); ++j)
            {
                dxdt[index_map_[*j]] += flux;
            }
        }
    }

protected:

    const std::vector<Species>& species_;
    const reaction_rule_container_type& reaction_rules_;
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

class ODESimulatorRatelow
    : public Simulator<NetworkModel, ODEWorld>
{
public:

    typedef Simulator<NetworkModel, ODEWorld> base_type;

public:

    ODESimulatorRatelow(
        boost::shared_ptr<NetworkModel> model,
        boost::shared_ptr<ODEWorld> world)
        : base_type(model, world), dt_(inf)
    {
        initialize();
    }

    ODESimulatorRatelow(boost::shared_ptr<ODEWorld> world)
        : base_type(world), dt_(inf)
    {
        initialize();
    }

    void initialize()
    {
        const std::vector<Species> species(model_->list_species());
        for (std::vector<Species>::const_iterator
                 i(species.begin()); i != species.end(); ++i)
        {
            if (!(*world_).has_species(*i))
            {
                (*world_).reserve_species(*i);
            }
        }
    }

    // SimulatorTraits

    Real t(void) const
    {
        return (*world_).t();
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

    Real dt_;
};

} // ode

} // ecell4

#endif /* __ECELL4_ODE_ODE_SIMULATOR_HPP */
