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
    typedef std::vector<state_type::size_type> index_container_type;

    typedef struct reaction
    {
        index_container_type reactants;
        index_container_type products;
        Real k;
    } reaction_type;

    typedef std::vector<reaction_type> reaction_container_type;

public:

    ODESystem(const reaction_container_type& reactions, const Real& volume)
        : reactions_(reactions), volume_(volume), vinv_(1.0 / volume)
    {
        ;
    }

    void operator()(const state_type& x, state_type& dxdt, const double& t)
    {
        for (state_type::iterator i(dxdt.begin()); i != dxdt.end(); ++i)
        {
            *i = 0.0;
        }

        for (reaction_container_type::const_iterator
            i(reactions_.begin()); i != reactions_.end(); ++i)
        {
            double flux((*i).k * volume_);

            for (index_container_type::const_iterator
                j((*i).reactants.begin()); j != (*i).reactants.end(); ++j)
            {
                flux *= x[*j] * vinv_;
            }

            for (index_container_type::const_iterator
                j((*i).reactants.begin()); j != (*i).reactants.end(); ++j)
            {
                dxdt[*j] -= flux;
            }

            for (index_container_type::const_iterator
                j((*i).products.begin()); j != (*i).products.end(); ++j)
            {
                dxdt[*j] += flux;
            }
        }
    }

protected:

    const reaction_container_type& reactions_;
    const Real volume_;
    const Real vinv_;
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
    : public Simulator<NetworkModel, ODEWorld>
{
public:

    typedef Simulator<NetworkModel, ODEWorld> base_type;

public:

    ODESimulator(
        boost::shared_ptr<NetworkModel> model,
        boost::shared_ptr<ODEWorld> world)
        : base_type(model, world), dt_(inf)
    {
        initialize();
    }

    ODESimulator(boost::shared_ptr<ODEWorld> world)
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

    ODESystem generate_system() const;

protected:

    Real dt_;
};

} // ode

} // ecell4

#endif /* __ECELL4_ODE_ODE_SIMULATOR_HPP */
