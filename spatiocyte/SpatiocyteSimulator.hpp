#ifndef __SPATIOCYTE_SIMULATOR_HPP
#define __SPATIOCYTE_SIMULATOR_HPP

#include <stdexcept>
#include <boost/shared_ptr.hpp>

#include <ecell4/core/Model.hpp>
#include <ecell4/core/Simulator.hpp>

#include "SpatiocyteWorld.hpp"


namespace ecell4
{

namespace spatiocyte
{

class SpatiocyteSimulator
    : public Simulator
{
public:

    SpatiocyteSimulator(boost::shared_ptr<Model> model, boost::shared_ptr<SpatiocyteWorld> world)
        : model_(model), world_(world), dt_(0), num_steps_(0)
    {
        ;
    }

    // SimulatorTraits

    Real t() const
    {
        return (*world_).t();
    }

    Real dt() const
    {
        return dt_;
    }

    Integer num_steps() const
    {
        return num_steps_;
    }

    void step();
    bool step(const Real& upto);

    // Optional members

    void set_t(const Real& t)
    {
        (*world_).set_t(t);
    }

    void set_dt(const Real& dt)
    {
        if (dt <= 0)
        {
            throw std::invalid_argument("The step size must be positive.");
        }
        dt_ = dt;
    }

    inline boost::shared_ptr<RandomNumberGenerator> rng()
    {
        return (*world_).rng();
    }

protected:

    boost::shared_ptr<Model> model_;
    boost::shared_ptr<SpatiocyteWorld> world_;

    /**
     * the protected internal state of SpatiocyteSimulator.
     * they are needed to be saved/loaded with Visitor pattern.
     */
    Real dt_;
    Integer num_steps_;
};

} // spatiocyte

} // ecell4

#endif /* __SPATIOCYTE_SIMULATOR_HPP */
