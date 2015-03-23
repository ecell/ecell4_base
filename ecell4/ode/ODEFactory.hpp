#ifndef __ECELL4_ODE_ODE_FACTORY_HPP
#define __ECELL4_ODE_ODE_FACTORY_HPP

#include <ecell4/core/SimulatorFactory.hpp>
#include <ecell4/core/extras.hpp>

#include "ODEWorld.hpp"
#include "ODESimulator.hpp"


namespace ecell4
{

namespace ode
{

class ODEFactory:
    public SimulatorFactory<ODEWorld, ODESimulator>
{
public:

    typedef SimulatorFactory<ODEWorld, ODESimulator> base_type;

public:

    ODEFactory()
        : base_type()
    {
        ; // do nothing
    }

    virtual ~ODEFactory()
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
        ODEWorld* w(extras::__generate_world_from_model(*this, m));

        for (Model::parameter_container_type::const_iterator
            i(m->parameters().begin()); i != m->parameters().end(); ++i)
        {
            const Species& sp(*i);
            if (sp.has_attribute("N"))
            {
                w->set_value(
                    Species(sp.serial()),
                    std::atof(sp.get_attribute("N").c_str())); //XXX: Real
            }
        }
        return w;
    }

    virtual ODESimulator* create_simulator(
        const boost::shared_ptr<Model>& model,
        const boost::shared_ptr<world_type>& world) const
    {
        return new ODESimulator(model, world);
        // return new ODESimulator(
        //     boost::dynamic_pointer_cast<NetworkModel>(model), world); //XXX:
    }

    virtual ODESimulator* create_simulator(
        const boost::shared_ptr<world_type>& world) const
    {
        return new ODESimulator(world);
    }

protected:
};

} // ode

} // ecell4

#endif /* __ECELL4_ODE_ODE_FACTORY_HPP */
