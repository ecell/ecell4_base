#ifndef __ODESIMULATOR_HPP
#define __ODESIMULATOR_HPP

#include <boost/shared_ptr.hpp>

#include <ecell4/core/exceptions.hpp>
#include <ecell4/core/types.hpp>
#include <ecell4/core/NetworkModel.hpp>
#include <ecell4/core/Simulator.hpp>

#include "ODEWorld.hpp"

namespace ecell4
{

namespace ode
{

class ODESystem
{
public:

    typedef std::vector<double> state_type;

    ODESystem()
    {
        ;
    }

    void operator()(state_type const& x, state_type& dxdt, double const& t)
    {
        ;
    }
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

    void operator()(ODESystem::state_type const&x, double t)
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
        : model_(model), world_(world), t_(0), num_steps_(0)
    {
        ;
    }

    void step(void)
    {
        throw NotImplemented("a step size must be specified.");
    }

    bool step(Real const& upto);

    Integer num_steps(void) const
    {
        return num_steps_;
    }

    Real t(void) const
    {
        return t_;
    }

protected:

    boost::shared_ptr<NetworkModel> model_;
    boost::shared_ptr<ODEWorld> world_;

    Real t_;
    Integer num_steps_;
};

} // ode

} // ecell4

#endif //__ODESIMULATOR_HPP
