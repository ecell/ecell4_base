#ifndef ECELL4_BD_BD_FACTORY_HPP
#define ECELL4_BD_BD_FACTORY_HPP

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

    BDFactory(const Integer3& matrix_sizes = default_matrix_sizes(), Real bd_dt_factor = default_bd_dt_factor())
        : base_type(), rng_(), matrix_sizes_(matrix_sizes), bd_dt_factor_(bd_dt_factor)
    {
        ; // do nothing
    }

    static inline const Integer3 default_matrix_sizes()
    {
        return Integer3(3, 3, 3);
    }

    static inline const Real default_bd_dt_factor()
    {
        return -1.0;
    }

    virtual ~BDFactory()
    {
        ; // do nothing
    }

    BDFactory& rng(const boost::shared_ptr<RandomNumberGenerator>& rng)
    {
        rng_ = rng;
        return (*this);
    }

    inline BDFactory* rng_ptr(const boost::shared_ptr<RandomNumberGenerator>& rng)
    {
        return &(this->rng(rng));  //XXX: == this
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
            return new BDWorld(edge_lengths, matrix_sizes_, rng_);
        }
        else
        {
            return new BDWorld(edge_lengths, matrix_sizes_);
        }
    }

    virtual BDWorld* create_world(const boost::shared_ptr<Model>& m) const
    {
        return extras::generate_world_from_model(*this, m);
    }

    virtual BDSimulator* create_simulator(
        const boost::shared_ptr<Model>& model,
        const boost::shared_ptr<world_type>& world) const
    {
        if (bd_dt_factor_ > 0)
        {
            return new BDSimulator(model, world, bd_dt_factor_);
        }
        else
        {
            return new BDSimulator(model, world);
        }
    }

    virtual BDSimulator* create_simulator(
        const boost::shared_ptr<world_type>& world) const
    {
        if (bd_dt_factor_ > 0)
        {
            return new BDSimulator(world, bd_dt_factor_);
        }
        else
        {
            return new BDSimulator(world);
        }
    }

protected:

    boost::shared_ptr<RandomNumberGenerator> rng_;
    Integer3 matrix_sizes_;
    Real bd_dt_factor_;
};

} // bd

} // ecell4

#endif /* ECELL4_BD_BD_FACTORY_HPP */
