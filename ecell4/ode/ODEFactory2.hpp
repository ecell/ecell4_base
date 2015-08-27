#ifndef __ECELL4_ODE_ODE_FACTORY2_HPP
#define __ECELL4_ODE_ODE_FACTORY2_HPP

#include <ecell4/core/SimulatorFactory.hpp>
#include <ecell4/core/extras.hpp>

#include "ODEWorld.hpp"
#include "ODESimulator2.hpp"


namespace ecell4
{

namespace ode
{

class ODEFactory2:
    public SimulatorFactory<ODEWorld, ODESimulator2>
{
public:

    typedef SimulatorFactory<ODEWorld, ODESimulator2> base_type;

public:

    ODEFactory2()
        : base_type()
    {
        ; // do nothing
    }

    virtual ~ODEFactory2()
    {
        ; // do nothing
    }

    virtual ODEWorld* create_world(const std::string filename) const
    {
        return new ODEWorld(filename);
    }

    virtual ODEWorld* create_world(
        const Real3& edge_lengths = Real3(1, 1, 1)) const
    {
        return new ODEWorld(edge_lengths);
    }

    virtual ODEWorld* create_world(const boost::shared_ptr<Model>& m) const
    {
        throw NotSupported("not supported.");
    }

    ODESimulator2* create_simulator(
        const boost::shared_ptr<Model>& model,
        const boost::shared_ptr<world_type>& world) const
    {
        throw NotSupported("not supported.");
    }

    virtual ODESimulator2* create_simulator(
        const boost::shared_ptr<world_type>& world) const
    {
        throw NotSupported("not supported.");
    }

    ODESimulator2* create_simulator(
        const boost::shared_ptr<NetworkModel>& model,
        const boost::shared_ptr<world_type>& world) const
    {
        return new ODESimulator2(model, world);
    }

    ODESimulator2* create_simulator(
        const boost::shared_ptr<ODENetworkModel>& model,
        const boost::shared_ptr<world_type>& world) const
    {
        return new ODESimulator2(model, world);
    }

protected:
};

} // ode

} // ecell4

#endif /* __ECELL4_ODE_ODE_FACTORY2_HPP */
