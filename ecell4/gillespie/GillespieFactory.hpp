#ifndef ECELL4_GILLESPIE_GILLESPIE_FACTORY_HPP
#define ECELL4_GILLESPIE_GILLESPIE_FACTORY_HPP

#include <ecell4/core/SimulatorFactory.hpp>
#include <ecell4/core/RandomNumberGenerator.hpp>

#include <ecell4/core/extras.hpp>
#include "GillespieWorld.hpp"
#include "GillespieSimulator.hpp"


namespace ecell4
{

namespace gillespie
{

class GillespieFactory:
    public SimulatorFactory<GillespieWorld, GillespieSimulator>
{
public:

    typedef SimulatorFactory<GillespieWorld, GillespieSimulator> base_type;
    typedef base_type::world_type world_type;
    typedef base_type::simulator_type simulator_type;
    typedef GillespieFactory this_type;

public:

    GillespieFactory()
        : base_type(), rng_()
    {
        ; // do nothing
    }

    virtual ~GillespieFactory()
    {
        ; // do nothing
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
            return new world_type(edge_lengths, rng_);
        }
        else
        {
            return new world_type(edge_lengths);
        }
    }

protected:

    boost::shared_ptr<RandomNumberGenerator> rng_;
};

} // gillespie

} // ecell4

#endif /* ECELL4_GILLESPIE_GILLESPIE_FACTORY_HPP */
