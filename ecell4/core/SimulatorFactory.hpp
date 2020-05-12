#ifndef ECELL4_SIMULATOR_FACTORY_HPP
#define ECELL4_SIMULATOR_FACTORY_HPP

#include "WorldInterface.hpp"
#include "Model.hpp"
#include "Simulator.hpp"
#include "extras.hpp"
#include "functions.hpp"


namespace ecell4
{

template <typename Tworld_, typename Tsim_>
class SimulatorFactory
{
public:

    typedef Tworld_ world_type;
    typedef Tsim_ simulator_type;

public:

    SimulatorFactory()
    {
        ; // do nothing
    }

    virtual ~SimulatorFactory()
    {
        ; // do nothing
    }

    world_type* world(const Real3& edge_lengths = ones()) const
    {
        return create_world(edge_lengths);
    }

    world_type* world(const Real volume) const
    {
        return world(ones() * cbrt(volume));
    }

    world_type* world(const std::string& filename) const
    {
        return new world_type(filename);
    }

    world_type* world(const std::shared_ptr<Model>& m) const
    {
        return extras::generate_world_from_model(*this, m);
    }

    simulator_type* simulator(
        const std::shared_ptr<world_type>& w, const std::shared_ptr<Model>& m) const
    {
        return create_simulator(w, m);
    }

    simulator_type* simulator(const std::shared_ptr<world_type>& w) const
    {
        if (std::shared_ptr<Model> bound_model = w->lock_model())
        {
            return create_simulator(w, bound_model);
        }
        else
        {
            throw std::invalid_argument("A world must be bound to a model.");
        }
    }

protected:

    virtual world_type* create_world(const Real3& edge_lengths) const
    {
        return new world_type(edge_lengths);
    }

    virtual simulator_type* create_simulator(
        const std::shared_ptr<world_type>& w, const std::shared_ptr<Model>& m) const
    {
        return new simulator_type(w, m);
    }
};

} // ecell4

#endif /* ECELL4_SIMULATOR_FACTORY_HPP */
