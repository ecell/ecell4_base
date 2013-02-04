#ifndef __ODESIMULATOR_HPP
#define __ODESIMULATOR_HPP

#include <boost/shared_ptr.hpp>

#include <ecell4/core/types.hpp>
#include <ecell4/core/NetworkModel.hpp>
#include <ecell4/core/Simulator.hpp>

#include "ODEWorld.hpp"

namespace ecell4
{

namespace ode
{

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
        ;
    }

    bool step(Real const& upto)
    {
        ;
    }

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
