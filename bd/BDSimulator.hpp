#ifndef __BD_SIMULATOR_HPP
#define __BD_SIMULATOR_HPP

#include <stdexcept>
#include <boost/shared_ptr.hpp>

#include <ecell4/core/RandomNumberGenerator.hpp>
#include <ecell4/core/Model.hpp>
#include <ecell4/core/Simulator.hpp>

#include "BDWorld.hpp"
#include "BDSimulatorState.hpp"
#include "BDPropagator.hpp"


namespace ecell4
{

namespace bd
{

class BDSimulator
    : public Simulator
{
public:

    BDSimulator(
        boost::shared_ptr<Model> model, boost::shared_ptr<BDWorld> world,
        RandomNumberGenerator& rng)
        : model_(model), world_(world), state_(new BDSimulatorState(rng))
    {
        ;
    }

    Real t() const
    {
        return (*world_).t();
    }

    void set_t(Real const& t)
    {
        (*world_).set_t(t);
    }

    Real dt() const
    {
        return (*state_).dt;
    }

    void set_dt(Real const& dt)
    {
        if (dt <= 0)
        {
            throw std::invalid_argument("The step size must be positive.");
        }
        (*state_).dt = dt;
    }

    Integer num_steps() const
    {
        return (*state_).num_steps;
    }

    void step();
    bool step(Real const& upto);

protected:

    boost::shared_ptr<Model> model_;
    boost::shared_ptr<BDWorld> world_;
    boost::shared_ptr<BDSimulatorState> state_;
};

} // bd

} // ecell4

#endif /* __BD_SIMULATOR_HPP */
