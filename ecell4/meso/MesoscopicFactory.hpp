#ifndef ECELL4_MESO_MESOSCOPIC_FACTORY_HPP
#define ECELL4_MESO_MESOSCOPIC_FACTORY_HPP

#include <ecell4/core/SimulatorFactory.hpp>
#include <ecell4/core/RandomNumberGenerator.hpp>
#include <ecell4/core/extras.hpp>

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

    MesoscopicFactory(const Integer3& matrix_sizes = default_matrix_sizes(), const Real subvolume_length = default_subvolume_length())
        : base_type(), rng_(), matrix_sizes_(matrix_sizes), subvolume_length_(subvolume_length)
    {
        ; // do nothing
    }

    static inline const Integer3 default_matrix_sizes()
    {
        return Integer3(1, 1, 1);
    }

    static inline const Real default_subvolume_length()
    {
        return 0.0;
    }

    virtual ~MesoscopicFactory()
    {
        ; // do nothing
    }

    MesoscopicFactory& rng(const boost::shared_ptr<RandomNumberGenerator>& rng)
    {
        rng_ = rng;
        return (*this);
    }

    inline MesoscopicFactory* rng_ptr(const boost::shared_ptr<RandomNumberGenerator>& rng)
    {
        return &(this->rng(rng));  //XXX: == this
    }

    virtual MesoscopicWorld* create_world(const std::string filename) const
    {
        return new MesoscopicWorld(filename);
    }

    virtual MesoscopicWorld* create_world(
        const Real3& edge_lengths = Real3(1, 1, 1)) const
    {
        if (rng_)
        {
            if (matrix_sizes_ != default_matrix_sizes())
            {
                return new MesoscopicWorld(edge_lengths, matrix_sizes_, rng_);
            }
            else if (subvolume_length_ != default_subvolume_length())
            {
                return new MesoscopicWorld(edge_lengths, subvolume_length_, rng_);
            }
            else
            {
                throw NotSupported(
                    "Either matrix_sizes or subvolume_length must be given.");
            }
        }
        if (matrix_sizes_ != default_matrix_sizes())
        {
            return new MesoscopicWorld(edge_lengths, matrix_sizes_);
        }
        else if (subvolume_length_ != default_subvolume_length())
        {
            return new MesoscopicWorld(edge_lengths, subvolume_length_);
        }
        return new MesoscopicWorld(edge_lengths);
    }

    virtual MesoscopicWorld* create_world(const boost::shared_ptr<Model>& m) const
    {
        return extras::generate_world_from_model(*this, m);
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

    boost::shared_ptr<RandomNumberGenerator> rng_;
    Integer3 matrix_sizes_;
    Real subvolume_length_;
};

} // meso

} // ecell4

#endif /* ECELL4_MESO_MESOSCOPIC_FACTORY_HPP */
