#ifndef ECELL4_SIMULATOR_FACTORY_HPP
#define ECELL4_SIMULATOR_FACTORY_HPP

#include "WorldInterface.hpp"
#include "Model.hpp"
#include "Simulator.hpp"
#include "extras.hpp"


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

    virtual world_type* world(const Real3& edge_lengths = ones()) const
    {
        return new world_type(edge_lengths);
    }

    virtual world_type* world(const std::string filename) const
    {
        return new world_type(filename);
    }

    virtual world_type* world(const boost::shared_ptr<Model>& m) const
    {
        return extras::generate_world_from_model(*this, m);
    }

    virtual simulator_type* simulator(
        const boost::shared_ptr<world_type>& w, const boost::shared_ptr<Model>& m) const
    {
        return new simulator_type(w, m);
    }

    virtual simulator_type* simulator(const boost::shared_ptr<world_type>& w) const
    {
        return new simulator_type(w);
    }

    world_type* create_world(const std::string filename) const
    {
        return this->world(filename);
    }

    world_type* create_world(const Real3& edge_lengths = ones()) const
    {
        return this->world(edge_lengths);
    }

    world_type* create_world(const boost::shared_ptr<Model>& m) const
    {
        return this->world(m);
    }

    simulator_type* create_simulator(
        const boost::shared_ptr<world_type>& w, const boost::shared_ptr<Model>& m) const
    {
        return this->simulator(w, m);
    }

    simulator_type* create_simulator(const boost::shared_ptr<world_type>& w) const
    {
        return this->simulator(w);
    }
};

} // ecell4

#endif /* ECELL4_SIMULATOR_FACTORY_HPP */
