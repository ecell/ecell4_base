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

    GillespieFactory& rng(const boost::shared_ptr<RandomNumberGenerator>& rng)
    {
        rng_ = rng;
        return (*this);
    }

    inline GillespieFactory* rng_ptr(const boost::shared_ptr<RandomNumberGenerator>& rng)
    {
        return &(this->rng(rng));  //XXX: == this
    }

    virtual GillespieWorld* create_world(const std::string filename) const
    {
        return new GillespieWorld(filename);
    }

    virtual GillespieWorld* create_world(
        const Real3& edge_lengths = Real3(1, 1, 1)) const
    {
        if (rng_)
        {
            return new GillespieWorld(edge_lengths, rng_);
        }
        else
        {
            return new GillespieWorld(edge_lengths);
        }
    }

    virtual GillespieWorld* create_world(const boost::shared_ptr<Model>& m) const
    {
        return extras::generate_world_from_model(*this, m);
    }

    virtual GillespieSimulator* create_simulator(
        const boost::shared_ptr<Model>& model,
        const boost::shared_ptr<world_type>& world) const
    {
        return new GillespieSimulator(model, world);
    }

    virtual GillespieSimulator* create_simulator(
        const boost::shared_ptr<world_type>& world) const
    {
        return new GillespieSimulator(world);
    }

protected:

    boost::shared_ptr<RandomNumberGenerator> rng_;
};

} // gillespie

} // ecell4

#endif /* ECELL4_GILLESPIE_GILLESPIE_FACTORY_HPP */
