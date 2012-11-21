#ifndef __BD_PROPAGATOR_HPP
#define __BD_PROPAGATOR_HPP

#include <ecell4/core/RandomNumberGenerator.hpp>
#include <ecell4/core/Model.hpp>

#include "functions3d.hpp"
#include "BDWorld.hpp"
#include "BDSimulatorState.hpp"


namespace ecell4
{

namespace bd
{

struct SpeciesInfo
{
    Real const radius;
    Real const D;
};

class BDPropagator
{
public:

    BDPropagator(Model& model, BDWorld& world, BDSimulatorState& state)
        : model_(model), world_(world), state_(state), max_retry_count_(1)
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

    class particle_finder
        : public std::unary_function<std::pair<ParticleID, Particle>, bool>
    {
    public:

        particle_finder(ParticleID const& pid)
            : pid_(pid)
        {
            ;
        }

        bool operator()(std::pair<ParticleID, Particle> pid_particle_pair)
        {
            return (pid_particle_pair.first == pid_);
        }

    protected:

        ParticleID pid_;
    };

    bool remove_particle(ParticleID const& pid);

    SpeciesInfo get_species_info(Species const& sp) const;

    inline Position3 draw_displacement(Particle const& particle)
    {
        return random_displacement_3d(rng(), dt(), particle.D());
    }

    inline Position3 draw_ipv(Real const& sigma, Real const& t, Real const& D)
    {
        return random_ipv_3d(rng(), sigma, t, D);
    }

protected:

    std::vector<std::pair<ParticleID, Particle> > queue_;
    Integer max_retry_count_;

    Model& model_;
    BDWorld& world_;
    BDSimulatorState& state_;
};

} // bd

} // ecell4

#endif /* __BD_PROPAGATOR_HPP */
