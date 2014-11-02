#ifndef __ECELL4_MESO_MESOSCOPIC_FACTORY_HPP
#define __ECELL4_MESO_MESOSCOPIC_FACTORY_HPP

#include <ecell4/core/SimulatorFactory.hpp>
#include <ecell4/core/RandomNumberGenerator.hpp>

#include "MesoscopicWorld.hpp"
#include "MesoscopicSimulator.hpp"


namespace ecell4
{

namespace meso
{

class MesoscopicFactory:
    public SimulatorFactory<MesoscopicWorld, MesoscopicSimulator>
{
public:

    typedef SimulatorFactory<MesoscopicWorld, MesoscopicSimulator> base_type;

public:

    MesoscopicFactory()
        : base_type(), matrix_sizes_(0, 0, 0), rng_()
    {
        ; // do nothing
    }

    MesoscopicFactory(const Integer3& matrix_sizes)
        : base_type(), matrix_sizes_(matrix_sizes), rng_()
    {
        ; // do nothing
    }

    MesoscopicFactory(const Integer3& matrix_sizes,
        const boost::shared_ptr<RandomNumberGenerator>& rng)
        : base_type(), matrix_sizes_(matrix_sizes), rng_(rng)
    {
        ; // do nothing
    }

    virtual ~MesoscopicFactory()
    {
        ; // do nothing
    }

    virtual MesoscopicWorld* create_world(const std::string filename) const
    {
        return new MesoscopicWorld(filename);
    }

    virtual MesoscopicWorld* create_world(
        const Position3& edge_lengths = Position3(1, 1, 1)) const
    {
        if (rng_)
        {
            return new MesoscopicWorld(edge_lengths, matrix_sizes_, rng_);
        }
        else if (matrix_sizes_[0] > 0 && matrix_sizes_[1] > 0 && matrix_sizes_[2] > 0)
        {
            return new MesoscopicWorld(edge_lengths, matrix_sizes_);
        }
        else
        {
            return new MesoscopicWorld(edge_lengths);
        }
    }

    virtual MesoscopicSimulator* create_simulator(
        const boost::shared_ptr<Model>& model,
        const boost::shared_ptr<world_type>& world) const
    {
        return new MesoscopicSimulator(model, world);
    }

    virtual MesoscopicSimulator* create_simulator(
        const boost::shared_ptr<world_type>& world) const
    {
        return new MesoscopicSimulator(world);
    }

protected:

    Integer3 matrix_sizes_;
    boost::shared_ptr<RandomNumberGenerator> rng_;
};

} // meso

} // ecell4

#endif /* __ECELL4_MESO_MESOSCOPIC_FACTORY_HPP */
