#ifndef ECELL4_LATTICE_LATTICE_FACTORY_HPP
#define ECELL4_LATTICE_LATTICE_FACTORY_HPP

#include <ecell4/core/RandomNumberGenerator.hpp>
#include <ecell4/core/SimulatorFactory.hpp>

#include "SpatiocyteSimulator.hpp"
#include "SpatiocyteWorld.hpp"
#include <ecell4/core/extras.hpp>

namespace ecell4
{

namespace spatiocyte
{

class SpatiocyteFactory
    : public SimulatorFactory<SpatiocyteWorld, SpatiocyteSimulator> {
public:
    typedef SimulatorFactory<SpatiocyteWorld, SpatiocyteSimulator> base_type;
    typedef base_type::world_type world_type;
    typedef base_type::simulator_type simulator_type;
    typedef SpatiocyteFactory this_type;

public:
    SpatiocyteFactory(const Real voxel_radius = default_voxel_radius())
        : base_type(), rng_(), voxel_radius_(voxel_radius)
    {
        ; // do nothing
    }

    virtual ~SpatiocyteFactory()
    {
        ; // do nothing
    }

    static inline const Real default_voxel_radius() { return 0.0; }

    this_type &rng(const boost::shared_ptr<RandomNumberGenerator> &rng)
    {
        rng_ = rng;
        return (*this);
    }

    inline this_type *
    rng_ptr(const boost::shared_ptr<RandomNumberGenerator> &rng)
    {
        return &(this->rng(rng)); // XXX: == this
    }

protected:
    virtual world_type *create_world(const Real3 &edge_lengths) const
    {
        if (rng_)
        {
            return new world_type(edge_lengths, voxel_radius_, rng_);
        }
        else if (voxel_radius_ > 0)
        {
            return new world_type(edge_lengths, voxel_radius_);
        }
        else
        {
            return new world_type(edge_lengths);
        }
    }

protected:
    boost::shared_ptr<RandomNumberGenerator> rng_;
    Real voxel_radius_;
};

} // namespace spatiocyte

} // namespace ecell4

#endif /* ECELL4_LATTICE_LATTICE_FACTORY_HPP */
