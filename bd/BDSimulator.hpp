#ifndef __BD_SIMULATOR_HPP
#define __BD_SIMULATOR_HPP

#include <stdexcept>
#include <boost/shared_ptr.hpp>

#include <ecell4/core/RandomNumberGenerator.hpp>
#include <ecell4/core/Model.hpp>
#include <ecell4/core/Simulator.hpp>

#include "BDWorld.hpp"


namespace ecell4
{

namespace bd
{

Real I_bd_3d(Real const& r01, Real const& dt, Real const& D);

struct BDSimulatorState
{
    BDSimulatorState(RandomNumberGenerator& r)
        : rng(r), num_steps(0), dt(0)
    {
        ;
    }

    Real dt;
    Integer num_steps;
    RandomNumberGenerator& rng;
};

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

    bool attempt_reaction(ParticleID const& pid, Particle const& particle);
    bool attempt_reaction(
        ParticleID const& pid1, Particle const& particle1,
        ParticleID const& pid2, Particle const& particle2);

    Position3 draw_displacement_3d(Particle const& particle);

    inline Position3 draw_displacement(Particle const& particle)
    {
        return draw_displacement_3d(particle);
    }

protected:

    boost::shared_ptr<Model> model_;
    boost::shared_ptr<BDWorld> world_;
    boost::shared_ptr<BDSimulatorState> state_;
};

} // bd

} // ecell4

#endif /* __BD_SIMULATOR_HPP */
