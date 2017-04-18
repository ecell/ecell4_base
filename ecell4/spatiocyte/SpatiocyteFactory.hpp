#ifndef ECELL4_LATTICE_LATTICE_FACTORY_HPP
#define ECELL4_LATTICE_LATTICE_FACTORY_HPP

#include <ecell4/core/SimulatorFactory.hpp>
#include <ecell4/core/RandomNumberGenerator.hpp>

#include <ecell4/core/extras.hpp>
#include "SpatiocyteWorld.hpp"
#include "SpatiocyteSimulator.hpp"


namespace ecell4
{

namespace spatiocyte
{

class SpatiocyteFactory:
    public SimulatorFactory<SpatiocyteWorld, SpatiocyteSimulator>
{
public:

    typedef SimulatorFactory<SpatiocyteWorld, SpatiocyteSimulator> base_type;

public:

    SpatiocyteFactory(const Real voxel_radius = default_voxel_radius())
        : base_type(), rng_(), voxel_radius_(voxel_radius)
    {
        ; // do nothing
    }

    static inline const Real default_voxel_radius()
    {
        return 0.0;
    }

    virtual ~SpatiocyteFactory()
    {
        ; // do nothing
    }

    SpatiocyteFactory& rng(const boost::shared_ptr<RandomNumberGenerator>& rng)
    {
        rng_ = rng;
        return (*this);
    }

    inline SpatiocyteFactory* rng_ptr(const boost::shared_ptr<RandomNumberGenerator>& rng)
    {
        return &(this->rng(rng));  //XXX: == this
    }

    virtual SpatiocyteWorld* create_world(const std::string filename) const
    {
        return new SpatiocyteWorld(filename);
    }

    virtual SpatiocyteWorld* create_world(
        const Real3& edge_lengths = Real3(1, 1, 1)) const
    {
        if (rng_)
        {
            return new SpatiocyteWorld(edge_lengths, voxel_radius_, rng_);
        }
        else if (voxel_radius_ > 0)
        {
            return new SpatiocyteWorld(edge_lengths, voxel_radius_);
        }
        else
        {
            return new SpatiocyteWorld(edge_lengths);
        }
    }

    virtual SpatiocyteWorld* create_world(const boost::shared_ptr<Model>& m) const
    {
        return extras::generate_world_from_model(*this, m);
    }

    virtual SpatiocyteSimulator* create_simulator(
        const boost::shared_ptr<Model>& model,
        const boost::shared_ptr<world_type>& world) const
    {
        return new SpatiocyteSimulator(model, world);
    }

    virtual SpatiocyteSimulator* create_simulator(
        const boost::shared_ptr<world_type>& world) const
    {
        return new SpatiocyteSimulator(world);
    }

protected:

    Real voxel_radius_;
    boost::shared_ptr<RandomNumberGenerator> rng_;
};

} // spatiocyte

} // ecell4

#endif /* ECELL4_LATTICE_LATTICE_FACTORY_HPP */
