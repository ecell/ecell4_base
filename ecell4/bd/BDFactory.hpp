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
    typedef base_type::world_type world_type;
    typedef base_type::simulator_type simulator_type;
    typedef BDFactory this_type;

public:

    BDFactory(const Integer3& matrix_sizes = default_matrix_sizes(), Real bd_dt_factor = default_bd_dt_factor())
        : base_type(), rng_(), matrix_sizes_(matrix_sizes), bd_dt_factor_(bd_dt_factor)
    {
        ; // do nothing
    }

    virtual ~BDFactory()
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

    this_type& rng(const boost::shared_ptr<RandomNumberGenerator>& rng)
    {
        rng_ = rng;
        return (*this);
    }

    inline this_type* rng_ptr(const boost::shared_ptr<RandomNumberGenerator>& rng)
    {
        return &(this->rng(rng));  //XXX: == this
    }

protected:

    virtual world_type* create_world(const Real3& edge_lengths) const
    {
        if (rng_)
        {
            return new world_type(edge_lengths, matrix_sizes_, rng_);
        }
        else
        {
            return new world_type(edge_lengths, matrix_sizes_);
        }
    }

    virtual simulator_type* create_simulator(
        const boost::shared_ptr<world_type>& w, const boost::shared_ptr<Model>& m) const
    {
        if (bd_dt_factor_ > 0)
        {
            return new simulator_type(w, m, bd_dt_factor_);
        }
        else
        {
            return new simulator_type(w, m);
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
