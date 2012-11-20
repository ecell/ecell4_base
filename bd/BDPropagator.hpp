#ifndef __BD_PROPAGATOR_HPP
#define __BD_PROPAGATOR_HPP

#include <ecell4/core/RandomNumberGenerator.hpp>
#include <ecell4/core/Model.hpp>

#include "BDWorld.hpp"
#include "BDSimulatorState.hpp"


namespace ecell4
{

namespace bd
{

Position3 random_displacement_3d(
    RandomNumberGenerator& rng, Real const& t, Real const& D);
Real I_bd_3d(Real const& r01, Real const& dt, Real const& D);

class BDPropagator
{
public:

    BDPropagator(Model& model, BDWorld& world, BDSimulatorState& state)
        : model_(model), world_(world), state_(state)
    {
        queue_ = world_.get_particles();
        shuffle(state_.rng, queue_);
    }

    bool operator()();

    inline Real dt() const
    {
        return state_.dt;
    }

    inline RandomNumberGenerator& rng()
    {
        return state_.rng;
    }

    bool attempt_reaction(ParticleID const& pid, Particle const& particle);
    bool attempt_reaction(
        ParticleID const& pid1, Particle const& particle1,
        ParticleID const& pid2, Particle const& particle2);

    inline Position3 draw_displacement(Particle const& particle)
    {
        return random_displacement_3d(rng(), dt(), particle.D());
    }

protected:

    std::vector<std::pair<ParticleID, Particle> > queue_;

    Model& model_;
    BDWorld& world_;
    BDSimulatorState& state_;
};

} // bd

} // ecell4

#endif /* __BD_PROPAGATOR_HPP */
