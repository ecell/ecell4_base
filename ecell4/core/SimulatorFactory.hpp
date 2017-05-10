#ifndef ECELL4_SIMULATOR_FACTORY_HPP
#define ECELL4_SIMULATOR_FACTORY_HPP

#include "Space.hpp"
#include "Model.hpp"
#include "Simulator.hpp"


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

    virtual Space* create_world(const std::string filename) const = 0;
    virtual Space* create_world(const Real3& edge_lengths) const = 0;

    virtual Space* create_world(const boost::shared_ptr<Model>& m) const
    {
        throw NotSupported("not supported yet");
    }

    virtual Simulator* create_simulator(
        const boost::shared_ptr<Model>& model,
        const boost::shared_ptr<world_type>& world) const = 0;
    virtual Simulator* create_simulator(
        const boost::shared_ptr<world_type>& world) const = 0;
};

} // ecell4

#endif /* ECELL4_SIMULATOR_FACTORY_HPP */
