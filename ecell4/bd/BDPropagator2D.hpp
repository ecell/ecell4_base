#ifndef __ECELL4_BD_BD_PROPAGATOR_HPP
#define __ECELL4_BD_BD_PROPAGATOR_HPP

#include <ecell4/core/RandomNumberGenerator.hpp>
#include <ecell4/core/Model.hpp>

#include "functions3d.hpp"
#include "BDPropagator.hpp"


namespace ecell4
{

namespace bd
{

class BDPropagator2D
{
public:

    typedef ReactionInfo reaction_info_type;

public:

    BDPropagator(
        Model& model, BDWorld& world, BDPolygon& poly, RandomNumberGenerator& rng, const Real& dt,
        std::vector<std::pair<ReactionRule, reaction_info_type> >& last_reactions)
        : model_(model), world_(world), poly_(poly), rng_(rng), dt_(dt),
        last_reactions_(last_reactions), max_retry_count_(1)
    {
        queue_ = poly_.list_particles();
        shuffle(rng_, queue_);
    }

    bool operator()();

    inline Real dt() const
    {
        return dt_;
    }

    inline RandomNumberGenerator& rng()
    {
        return rng_;
    }

    bool attempt_reaction(const ParticleID& pid, const Particle& particle);
    bool attempt_reaction(
        const ParticleID& pid1, const Particle& particle1,
        const ParticleID& pid2, const Particle& particle2);

    class particle_finder
        : public std::unary_function<std::pair<ParticleID, Particle>, bool>
    {
    public:

        particle_finder(const ParticleID& pid)
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

    void remove_particle(const ParticleID& pid);

    inline Real3 draw_displacement(const Particle& particle)
    {// XXX!
        return random_displacement_2d(rng(), dt(), particle.D());
    }

    inline Real3 draw_ipv(const Real& sigma, const Real& t, const Real& D)
    {// XXX!
        return random_ipv_3d(rng(), sigma, t, D);
    }

protected:

    Model& model_;
    BDWorld& world_;
    BDPolygon& poly_;
    RandomNumberGenerator& rng_;
    Real dt_;
    std::vector<std::pair<ReactionRule, reaction_info_type> >& last_reactions_;
    Integer max_retry_count_;

    BDWorld::particle_container_type queue_;
};

} // bd

} // ecell4

#endif /* __ECELL4_BD_BD_PROPAGATOR_HPP */
