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
    typedef base_type::world_type world_type;
    typedef base_type::simulator_type simulator_type;
    typedef MesoscopicFactory this_type;

public:

    MesoscopicFactory(
        const Integer3& matrix_sizes = default_matrix_sizes(),
        const Real subvolume_length = default_subvolume_length())
        : base_type(), rng_(), matrix_sizes_(matrix_sizes), subvolume_length_(subvolume_length)
    {
        ; // do nothing
    }

    virtual ~MesoscopicFactory()
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

    this_type& rng(const boost::shared_ptr<RandomNumberGenerator>& rng)
    {
        rng_ = rng;
        return (*this);
    }

    inline this_type* rng_ptr(const boost::shared_ptr<RandomNumberGenerator>& rng)
    {
        return &(this->rng(rng));  //XXX: == this
    }

    virtual world_type* world(
        const Real3& edge_lengths = Real3(1, 1, 1)) const
    {
        if (rng_)
        {
            if (matrix_sizes_ != default_matrix_sizes())
            {
                return new world_type(edge_lengths, matrix_sizes_, rng_);
            }
            else if (subvolume_length_ != default_subvolume_length())
            {
                return new world_type(edge_lengths, subvolume_length_, rng_);
            }
            else
            {
                throw NotSupported(
                    "Either matrix_sizes or subvolume_length must be given.");
            }
        }
        if (matrix_sizes_ != default_matrix_sizes())
        {
            return new world_type(edge_lengths, matrix_sizes_);
        }
        else if (subvolume_length_ != default_subvolume_length())
        {
            return new world_type(edge_lengths, subvolume_length_);
        }
        return new world_type(edge_lengths);
    }

protected:

    boost::shared_ptr<RandomNumberGenerator> rng_;
    Integer3 matrix_sizes_;
    Real subvolume_length_;
};

} // meso

} // ecell4

#endif /* ECELL4_MESO_MESOSCOPIC_FACTORY_HPP */
