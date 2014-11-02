#ifndef __ECELL4_BD_BD_FACTORY_HPP
#define __ECELL4_BD_BD_FACTORY_HPP

#include <ecell4/core/SimulatorFactory.hpp>
#include <ecell4/core/RandomNumberGenerator.hpp>

#include "BDWorld.hpp"
#include "BDSimulator.hpp"


namespace ecell4
{

namespace bd
{

class BDFactory:
    public SimulatorFactory<BDWorld, BDSimulator>
{
public:

    typedef SimulatorFactory<BDWorld, BDSimulator> base_type;

public:

    BDFactory()
        : base_type(), rng_()
    {
        ; // do nothing
    }

    BDFactory(const boost::shared_ptr<RandomNumberGenerator>& rng)
        : base_type(), rng_(rng)
    {
        ; // do nothing
    }

    virtual ~BDFactory()
    {
        ; // do nothing
    }

    virtual BDWorld* create_world(const std::string filename) const
    {
        return new BDWorld(filename);
    }

    virtual BDWorld* create_world(
        const Real3& edge_lengths = Real3(1, 1, 1)) const
    {
        if (rng_)
        {
            return new BDWorld(edge_lengths, rng_);
        }
        else
        {
            return new BDWorld(edge_lengths);
        }
    }

    virtual BDSimulator* create_simulator(
        const boost::shared_ptr<Model>& model,
        const boost::shared_ptr<world_type>& world) const
    {
        return new BDSimulator(model, world);
    }

    virtual BDSimulator* create_simulator(
        const boost::shared_ptr<world_type>& world) const
    {
        return new BDSimulator(world);
    }

protected:

    boost::shared_ptr<RandomNumberGenerator> rng_;
};

} // bd

} // ecell4

#endif /* __ECELL4_BD_BD_FACTORY_HPP */
