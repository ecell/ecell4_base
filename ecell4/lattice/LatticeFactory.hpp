#ifndef __ECELL4_LATTICE_LATTICE_FACTORY_HPP
#define __ECELL4_LATTICE_LATTICE_FACTORY_HPP

#include <ecell4/core/SimulatorFactory.hpp>
#include <ecell4/core/RandomNumberGenerator.hpp>

#include "LatticeWorld.hpp"
#include "LatticeSimulator.hpp"


namespace ecell4
{

namespace lattice
{

class LatticeFactory:
    public SimulatorFactory<LatticeWorld, LatticeSimulator>
{
public:

    typedef SimulatorFactory<LatticeWorld, LatticeSimulator> base_type;

public:

    LatticeFactory()
        : base_type(), voxel_radius_(0.0), rng_()
    {
        ; // do nothing
    }

    LatticeFactory(const Real voxel_radius)
        : base_type(), voxel_radius_(voxel_radius), rng_()
    {
        ; // do nothing
    }

    LatticeFactory(const Real voxel_radius,
        const boost::shared_ptr<RandomNumberGenerator>& rng)
        : base_type(), voxel_radius_(voxel_radius), rng_(rng)
    {
        ; // do nothing
    }

    virtual ~LatticeFactory()
    {
        ; // do nothing
    }

    virtual LatticeWorld* create_world(const std::string filename) const
    {
        return new LatticeWorld(filename);
    }

    virtual LatticeWorld* create_world(
        const Position3& edge_lengths = Position3(1, 1, 1)) const
    {
        if (rng_)
        {
            return new LatticeWorld(edge_lengths, voxel_radius_, rng_);
        }
        else if (voxel_radius_ > 0)
        {
            return new LatticeWorld(edge_lengths, voxel_radius_);
        }
        else
        {
            return new LatticeWorld(edge_lengths);
        }
    }

    virtual LatticeSimulator* create_simulator(
        const boost::shared_ptr<Model>& model,
        const boost::shared_ptr<world_type>& world) const
    {
        return new LatticeSimulator(model, world);
    }

    virtual LatticeSimulator* create_simulator(
        const boost::shared_ptr<world_type>& world) const
    {
        return new LatticeSimulator(world);
    }

protected:

    Real voxel_radius_;
    boost::shared_ptr<RandomNumberGenerator> rng_;
};

} // lattice

} // ecell4

#endif /* __ECELL4_LATTICE_LATTICE_FACTORY_HPP */
